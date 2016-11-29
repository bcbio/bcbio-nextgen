"""Cohort based copy number calling in gene regions using Seq2C.

Seq2C calls across multiple samples without explicit background samples,
using gene regions as segments.

This requires coverage calculation in each sample and gene, followed by global
calling across all samples.
"""
import os
import subprocess
from collections import defaultdict
import pybedtools as bt

from bcbio import utils
from bcbio.variation.vcfutils import get_paired_phenotype
from bcbio.pipeline import datadict as dd, config_utils
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.structural import annotate, regions
from bcbio.log import logger
from bcbio.variation.coverage import regions_coverage
from bcbio.variation import bedutils


def precall(items):
    """Perform initial pre-calling steps -- coverage calcuation by sample.

    Use sambamba to call average region coverage in regions, and convert into a correct format.
    """
    items = [utils.to_single_data(x) for x in items]
    assert len(items) == 1, "Expect one item to Seq2C coverage calculation"
    data = utils.to_single_data(items)
    # sv_bed could specify a smaller region than variant coverage, so avoid
    # this sanity check
    # assert dd.get_coverage_interval(data) != "genome", "Seq2C only for amplicon and exome sequencing"

    assert "seq2c_bed_ready" in data["config"]["algorithm"], "Error: svregions or variant_regions BED file required for Seq2C"

    bed_file = data["config"]["algorithm"]["seq2c_bed_ready"]
    bam_file = dd.get_align_bam(data)
    sample_name = dd.get_sample_name(data)

    work_dir = _sv_workdir(data)
    cov_file = _calculate_coverage(data, work_dir, bed_file, bam_file, sample_name)

    if "sv" not in data:
        data["sv"] = []
    data["sv"].append({"variantcaller": "seq2c",
                       "coverage": cov_file})
    return [data]

def run(items):
    """Normalization and log2 ratio calculation plus CNV calling for full cohort.

    - Combine coverage of each region for each sample
    - Prepare read counts for each sample
    - Normalize coverages in cohort by gene and sample, and calculate log2 ratios
    - Call amplifications and deletions
    """
    items = [utils.to_single_data(x) for x in items]
    work_dir = _sv_workdir(items[0])

    coverage_file = _combine_coverages(items, work_dir)
    read_mapping_file = _calculate_mapping_reads(items, work_dir)

    normal_names = [dd.get_sample_name(x) for x in items if get_paired_phenotype(x) == "normal"]
    seq2c_calls_file = _call_cnv(items, work_dir, read_mapping_file, coverage_file, normal_names)
    _split_cnv(items, seq2c_calls_file)

    return items

def prep_seq2c_bed(data):
    """Selecting the bed file, cleaning, and properly annotating for Seq2C
    """
    bed_file = regions.get_sv_bed(data)
    if bed_file:
        bed_file = bedutils.clean_file(bed_file, data, prefix="svregions-")
    else:
        bed_file = bedutils.clean_file(dd.get_variant_regions(data), data)
    if not bed_file:
        return None

    col_num = bt.BedTool(bed_file).field_count()
    if col_num < 4:
        annotated_file = annotate.add_genes(bed_file, data, max_distance=0)
        if annotated_file == bed_file:
            raise ValueError("BED file for Seq2C must be annotated with gene names, "
                             "however the input BED is 3-columns and we have no transcript "
                             "data to annotate with " + bed_file)
        annotated_file = annotate.gene_one_per_line(annotated_file, data)
    else:
        annotated_file = bed_file

    ready_file = "%s-seq2cclean.bed" % (utils.splitext_plus(annotated_file)[0])
    if not utils.file_uptodate(ready_file, annotated_file):
        bed = bt.BedTool(annotated_file)
        if col_num > 4 and col_num != 8:
            bed = bed.cut(range(4))
        bed = bed.filter(lambda x: x.name not in ["", ".", "-"])
        with file_transaction(data, ready_file) as tx_out_file:
            bed.saveas(tx_out_file)
        logger.debug("Saved Seq2C clean annotated ready input BED into " + ready_file)

    return ready_file

def _call_cnv(items, work_dir, read_mapping_file, coverage_file, control_sample_names):
    output_fpath = os.path.join(work_dir, "calls_combined.tsv")
    cov2lr = "cov2lr.pl"
    lr2gene = "lr2gene.pl"
    control_opt = ""
    lr2gene_opt = ""
    if control_sample_names:
        control_opt = "-c " + ":".join(control_sample_names)
        lr2gene_opt = "-c"

    if not utils.file_exists(output_fpath):
        with file_transaction(items[0], output_fpath) as tx_out_file:
            export = utils.local_path_export()
            cmd = ("{export} {cov2lr} -a {control_opt} {read_mapping_file} {coverage_file} | " +
                   "{lr2gene} {lr2gene_opt} > {output_fpath}")
            do.run(cmd.format(**locals()), "Seq2C CNV calling")
    return output_fpath

def _split_cnv(items, calls_fpath):
    for item in items:
        if get_paired_phenotype(item) == "normal":
            continue

        sample_name = dd.get_sample_name(item)
        work_dir = _sv_workdir(item)
        out_fname = os.path.join(work_dir, sample_name + '-calls.tsv')
        if not utils.file_exists(out_fname):
            with file_transaction(item, out_fname) as tx:
                with open(tx, "w") as out, open(calls_fpath) as inp:
                    out.write(next(inp))
                    for l in inp:
                        if l.split("\t")[0] == sample_name:
                            out.write(l)
        item["sv"][0]["calls"] = out_fname

def _calculate_coverage(data, work_dir, bed_file, bam_file, sample_name):
    sambamba_depth_file = regions_coverage(data, bed_file, bam_file, "sv_regions")

    out_file = os.path.join(work_dir, sample_name + '-coverage.tsv')
    if not utils.file_exists(out_file):
        logger.debug('Converting sambamba depth output to cov2lr.pl input in ' + sample_name)
        with file_transaction(data, out_file) as tx_out_file:
            _sambabma_depth_to_seq2cov(sambamba_depth_file, tx_out_file, sample_name)
    logger.debug("Saved to " + out_file)
    return out_file

def _sambabma_depth_to_seq2cov(input_fpath, output_fpath, sample_name):
    """Args:
        input_fpath: output of "sambabma depth region":
            # chrom chromStart  chromEnd  F3       readCount  meanCoverage  sampleName
            chr20   68345       68413     DEFB125  56         32.5          chr20_tumor_1
            chr20   76640       77301     DEFB125  279        36.9213       chr20_tumor_1

        output_fpath: path to write results - input for Seq2C's cov2lr.pl, e.g.:
            seq2cov:
            chr20_tumor_1   DEFB125   chr20   68346   68413   Amplicon    68   28.0
            chr20_tumor_1   DEFB125   chr20   76641   77301   Amplicon    661  24.0
            chr20_tumor_1   DEFB125   chr20   68346   77301   Whole-Gene  729  24.3731138546

        sample_name:
            sample name (e.g. chr20_tumor_1)
    """
    # First round: collecting gene ends
    gene_end_by_gene = defaultdict(lambda: -1)
    with open(input_fpath) as f:
        ave_depth_col = next(f).split('\t').index('meanCoverage')
        for l in f:
            if l.startswith('#'): continue
            fs = l.replace('\n', '').split('\t')
            if any(fs[i] == '.' for i in [0, 1, 2, 3, ave_depth_col]): continue
            end = int(fs[2])
            gene_name = fs[3]
            gene_end_by_gene[gene_name] = max(gene_end_by_gene[gene_name], end)

    # Second round: calculating gene level coverage, and writing file for Seq2C
    total_cov_by_gene = dict()
    gene_start_by_gene = dict()
    total_size_by_gene = dict()

    with open(input_fpath) as f, open(output_fpath, 'w') as out:
        for l in f:
            if l.startswith('#'): continue
            fs = l.replace('\n', '').split('\t')
            if any(fs[i] == '.' for i in [0, 1, 2, 3, ave_depth_col]): continue

            chrom = fs[0]
            start = int(fs[1])
            end = int(fs[2])
            gene_name = fs[3]
            ave_depth = float(fs[ave_depth_col])

            if gene_name not in gene_start_by_gene:
                gene_start_by_gene[gene_name] = start
                total_cov_by_gene[gene_name] = 0
                total_size_by_gene[gene_name] = 0
            else:
                gene_start_by_gene[gene_name] = min(start, gene_start_by_gene[gene_name])
            total_cov_by_gene[gene_name] += ave_depth * (end - start)
            total_size_by_gene[gene_name] += end - start

            fs = [sample_name, gene_name, chrom, str(start + 1), str(end), 'Amplicon', str(end - start), str(ave_depth)]
            out.write('\t'.join(fs) + '\n')

            if end >= gene_end_by_gene[gene_name]:
                assert end == gene_end_by_gene[gene_name], (end, gene_end_by_gene[gene_name])
                start = gene_start_by_gene[gene_name]
                ave_depth = total_cov_by_gene[gene_name] / total_size_by_gene[gene_name]
                size = total_size_by_gene[gene_name]
                fs = [sample_name, gene_name, chrom, str(start + 1), str(end), 'Whole-Gene', str(size), str(ave_depth)]
                out.write('\t'.join(fs) + '\n')
    return output_fpath

def _combine_coverages(items, work_dir):
    """Combine coverage cnns calculated for individual inputs into single file.
    """
    out_file = os.path.join(work_dir, "sample_coverages.txt")
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            with open(tx_out_file, 'w') as out_f:
                for data in items:
                    svouts = [x for x in data["sv"] if x["variantcaller"] == "seq2c"]
                    assert len(svouts) == 1
                    cov_file = svouts[0]["coverage"]
                    with open(cov_file) as cov_f:
                        out_f.write(cov_f.read())
    return out_file

def _calculate_mapping_reads(items, work_dir):
    """Calculate read counts from samtools idxstats for each sample.
    """
    out_file = os.path.join(work_dir, "mapping_reads.txt")
    if not utils.file_exists(out_file):
        lines = []
        for data in items:
            count = 0
            for line in subprocess.check_output([
                "samtools", "idxstats", dd.get_align_bam(data)]).split("\n"):
                if line.strip():
                    count += int(line.split("\t")[2])
            lines.append("%s\t%s" % (dd.get_sample_name(data), count))
        with file_transaction(items[0], out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                out_handle.write("\n".join(lines))
    return out_file

# def _count_mapped_reads(data, work_dir, bed_file, bam_file):
#     """Calculate read counts from samtools idxstats for each sample.
#     """
#     sambamba = config_utils.get_program("sambamba", data["config"])
#     num_cores = dd.get_cores(data)
#     out_file = os.path.join(work_dir, "mapped_reads_count.txt")
#     if not utils.file_exists(out_file):
#         logger.debug('Counting mapped reads in ' + dd.get_sample_name(data))
#         with file_transaction(data, out_file) as tx_out_file:
#             cmd = ("{sambamba} view -c -t {num_cores} "
#                    "-F \"not duplicate and not failed_quality_control\" "
#                    "-L {bed_file} {bam_file} > {tx_out_file}")
#             do.run(cmd.format(**locals()), "Counting mapped reads")
#     return out_file
#
# def _combine_read_counts(items, work_dir):
#     out_file = os.path.join(work_dir, "mapped_reads_counts.txt")
#     if not utils.file_exists(out_file):
#         with file_transaction(items[0], out_file) as tx_out_file:
#             with open(tx_out_file, "w") as out_handle:
#                 for data in items:
#                     with open(data["sv"][0]["read_count"]) as f:
#                         out_handle.write(dd.get_sample_name(data) + "\t" + f.read())
#     return out_file

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "seq2c"))
