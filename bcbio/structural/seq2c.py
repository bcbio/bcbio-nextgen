"""Cohort based copy number calling in gene regions using Seq2C.

Seq2C calls across multiple samples without explicit background samples,
using gene regions as segments.

This requires coverage calculation in each sample and gene, followed by global
calling across all samples.
"""
import os
import shutil
import subprocess
from collections import defaultdict

import pybedtools as bt
import toolz as tz

from bcbio import utils
from bcbio.variation.vcfutils import get_paired_phenotype
from bcbio.pipeline import datadict as dd, config_utils
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.structural import annotate, regions
from bcbio.log import logger
from bcbio.variation.coverage import regions_coverage
from bcbio.variation import bedutils, population, vcfutils


def precall(data):
    """Perform initial pre-calling steps -- coverage calcuation by sample.

    Use mosdepth to call average region coverage in regions, and convert
    into seq2c format.
    """
    data = utils.to_single_data(data)
    bed_file = tz.get_in(["config", "algorithm", "seq2c_bed_ready"], data)
    if not bed_file:
        raise ValueError("Error: svregions or variant_regions BED file required for Seq2C")
    sample_name = dd.get_sample_name(data)
    work_dir = _sv_workdir(data)
    return _calculate_coverage(data, work_dir, bed_file, sample_name)

def run(items):
    """Normalization and log2 ratio calculation plus CNV calling for full cohort.

    - Combine coverage of each region for each sample
    - Prepare read counts for each sample
    - Normalize coverages in cohort by gene and sample, and calculate log2 ratios
    - Call amplifications and deletions
    """
    items = [utils.to_single_data(x) for x in items]
    work_dir = _sv_workdir(items[0])

    input_backs = list(set(filter(lambda x: x is not None,
                                  [dd.get_background_cnv_reference(d, "seq2c") for d in items])))
    coverage_file = _combine_coverages(items, work_dir, input_backs)
    read_mapping_file = _calculate_mapping_reads(items, work_dir, input_backs)
    normal_names = []
    if input_backs:
        with open(input_backs[0]) as in_handle:
            for line in in_handle:
                if len(line.split()) == 2:
                    normal_names.append(line.split()[0])
    normal_names += [dd.get_sample_name(x) for x in items if population.get_affected_status(x) == 1]
    seq2c_calls_file = _call_cnv(items, work_dir, read_mapping_file, coverage_file, normal_names)
    items = _split_cnv(items, seq2c_calls_file, read_mapping_file, coverage_file)
    return items

def prep_seq2c_bed(data):
    """Selecting the bed file, cleaning, and properly annotating for Seq2C
    """
    if dd.get_background_cnv_reference(data, "seq2c"):
        bed_file = _background_to_bed(dd.get_background_cnv_reference(data, "seq2c"), data)
    else:
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

def _background_to_bed(back_file, data):
    """Convert a seq2c background file with calls into BED regions for coverage.

    seq2c background files are a concatenation of mapping and sample_coverages from
    potentially multiple samples. We use the coverage information from the first
    sample to translate into BED.
    """
    out_file = os.path.join(utils.safe_makedir(os.path.join(dd.get_work_dir(data), "bedprep")),
                            "%s-regions.bed" % utils.splitext_plus(os.path.basename(back_file))[0])
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(back_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    sample = in_handle.readline().split("\t")[0]
                    for line in in_handle:
                        if line.startswith(sample) and len(line.split()) >= 5:
                            _, gene, chrom, start, end = line.split()[:5]
                            out_handle.write("%s\n" % ("\t".join([chrom, str(int(start) - 1), end, gene])))
    return out_file

def _get_seq2c_options(data):
    """Get adjustable, through resources, or default options for seq2c.
    """
    cov2lr_possible_opts = ["-F"]
    defaults = {}
    ropts = config_utils.get_resources("seq2c", data["config"]).get("options", [])
    assert len(ropts) % 2 == 0, "Expect even number of options for seq2c" % ropts
    defaults.update(dict(tz.partition(2, ropts)))
    cov2lr_out, lr2gene_out = [], []
    for k, v in defaults.items():
        if k in cov2lr_possible_opts:
            cov2lr_out += [str(k), str(v)]
        else:
            lr2gene_out += [str(k), str(v)]
    return cov2lr_out, lr2gene_out

def _call_cnv(items, work_dir, read_mapping_file, coverage_file, control_sample_names):
    output_fpath = os.path.join(work_dir, "calls_combined.tsv")
    cov2lr = "cov2lr.pl"
    lr2gene = "lr2gene.pl"
    cov2lr_opts, lr2gene_opts = _get_seq2c_options(items[0])
    if control_sample_names:
        cov2lr_opts += ["-c", ":".join(control_sample_names)]
        if "-c" not in lr2gene_opts:
            lr2gene_opts += ["-c"]
    cov2lr_opt = " ".join(cov2lr_opts)
    lr2gene_opt = " ".join(lr2gene_opts)

    if not utils.file_exists(output_fpath):
        with file_transaction(items[0], output_fpath) as tx_out_file:
            with utils.chdir(work_dir):
                export = utils.local_path_export()
                cmd = ("{export} {cov2lr} -a {cov2lr_opt} {read_mapping_file} {coverage_file} | " +
                    "{lr2gene} {lr2gene_opt} > {tx_out_file}")
                do.run(cmd.format(**locals()), "Seq2C CNV calling")
    return output_fpath

def _split_cnv(items, calls_fpath, read_mapping_file, coverage_file):
    out_items = []
    for item in items:
        cur_sv = {"variantcaller": "seq2c", "coverage": tz.get_in(["depth", "bins", "seq2c"], item)}
        if not get_paired_phenotype(item) == "normal":
            sample_name = dd.get_sample_name(item)
            work_dir = _sv_workdir(item)
            out_fname = os.path.join(work_dir, sample_name + '-calls.tsv')
            gender_fname = os.path.join(work_dir, 'gender_predicted.txt')
            out_gender_name = os.path.join(work_dir, sample_name + '-gender_predicted.txt')
            if not utils.file_exists(out_fname):
                with file_transaction(item, out_fname) as tx:
                    with open(tx, "w") as out, open(calls_fpath) as inp:
                        out.write(next(inp))
                        for l in inp:
                            if l.split("\t")[0] == sample_name:
                                out.write(l)
            cur_sv.update({"calls": out_fname, "vrn_file": to_vcf(out_fname, item),
                           "read_mapping": read_mapping_file, "calls_all": calls_fpath,
                           "coverage_all": coverage_file})
            if utils.file_exists(gender_fname):
                if not utils.file_exists(out_gender_name):
                    shutil.copyfile(gender_fname, out_gender_name)
                cur_sv.update({"gender_predicted": out_gender_name})
        if "sv" not in item:
            item["sv"] = []
        assert "seq2c" not in [x["variantcaller"] for x in item["sv"]], \
            "Do not expect existing seq2c variant output: %s" % (dd.get_sample_name(item))
        item["sv"].append(cur_sv)
        out_items.append(item)
    return out_items

VCF_HEADER = """##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=FOLD_CHANGE_LOG,Number=1,Type=Float,Description="Log fold change">
##INFO=<ID=PROBES,Number=1,Type=Integer,Description="Number of exons in change">
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene identified">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
"""

def to_vcf(in_tsv, data):
    """Convert seq2c output file into BED output.
    """
    call_convert = {"Amp": "DUP", "Del": "DEL"}
    out_file = "%s.vcf" % utils.splitext_plus(in_tsv)[0]
    if not utils.file_uptodate(out_file, in_tsv):
        with file_transaction(data, out_file) as tx_out_file:
            with open(in_tsv) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    out_handle.write(VCF_HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n"
                                     % (dd.get_sample_name(data)))
                    header = in_handle.readline().split("\t")
                    for cur in (dict(zip(header, l.split("\t"))) for l in in_handle):
                        if cur["Amp_Del"] in call_convert:
                            svtype = call_convert[cur["Amp_Del"]]
                            info = "SVTYPE=%s;END=%s;SVLEN=%s;FOLD_CHANGE_LOG=%s;PROBES=%s;GENE=%s" % (
                                svtype, cur["End"], int(cur["End"]) - int(cur["Start"]),
                                cur["Log2ratio"], cur["Ab_Seg"], cur["Gene"])
                            out_handle.write("\t".join([cur["Chr"], cur["Start"], ".", "N", "<%s>" % (svtype),
                                                        ".", ".", info, "GT", "1/1"]) + "\n")
    return vcfutils.sort_by_ref(out_file, data)

def _calculate_coverage(data, work_dir, bed_file, sample_name):
    depth_file = regions_coverage(bed_file, "seq2c_bed_ready", data)

    out_file = os.path.join(work_dir, sample_name + '-coverage.tsv')
    if not utils.file_exists(out_file):
        logger.debug('Converting depth output to cov2lr.pl input in ' + sample_name)
        with file_transaction(data, out_file) as tx_out_file:
            _depth_to_seq2cov(depth_file, tx_out_file, sample_name)
    logger.debug("Saved to " + out_file)
    return out_file

def _depth_to_seq2cov(input_fpath, output_fpath, sample_name):
    """Args:
        input_fpath: output of "mosdepth":
            chr22           14250   15500   name3   5.54
            chrM            100     1000    name1   916.08

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
    with utils.open_gzipsafe(input_fpath) as f:
        for xs in (l.rstrip().split() for l in f if not l.startswith("#")):
            xs = [x for x in xs if x.strip()]
            if any(x == "." for x in xs): continue
            end = int(xs[2])
            gene_name = xs[3]
            gene_end_by_gene[gene_name] = max(gene_end_by_gene[gene_name], end)

    # Second round: calculating gene level coverage, and writing file for Seq2C
    total_cov_by_gene = dict()
    gene_start_by_gene = dict()
    total_size_by_gene = dict()

    with utils.open_gzipsafe(input_fpath) as f, open(output_fpath, 'w') as out:
        for xs in (l.rstrip().split() for l in f if not l.startswith("#")):
            xs = [x for x in xs if x.strip()]
            if any(x == "." for x in xs): continue
            chrom, start, end, gene_name = xs[:4]
            start, end = int(start), int(end)
            ave_depth = float(xs[-1])

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

def _combine_coverages(items, work_dir, input_backs=None):
    """Combine coverage cnns calculated for individual inputs into single file.

    Optionally moves over pre-calculated coverage samples from a background file.
    """
    out_file = os.path.join(work_dir, "sample_coverages.txt")
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            with open(tx_out_file, 'w') as out_f:
                for data in items:
                    cov_file = tz.get_in(["depth", "bins", "seq2c"], data)
                    with open(cov_file) as cov_f:
                        out_f.write(cov_f.read())
                if input_backs:
                    for input_back in input_backs:
                        with open(input_back) as in_handle:
                            for line in in_handle:
                                if len(line.split()) >= 4:
                                    out_f.write(line)
    return out_file

def _calculate_mapping_reads(items, work_dir, input_backs=None):
    """Calculate read counts from samtools idxstats for each sample.

    Optionally moves over pre-calculated mapping counts from a background file.
    """
    out_file = os.path.join(work_dir, "mapping_reads.txt")
    if not utils.file_exists(out_file):
        lines = []
        for data in items:
            count = 0
            for line in subprocess.check_output([
                "samtools", "idxstats", dd.get_align_bam(data)]).decode().split("\n"):
                if line.strip():
                    count += int(line.split("\t")[2])
            lines.append("%s\t%s" % (dd.get_sample_name(data), count))
        with file_transaction(items[0], out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                out_handle.write("\n".join(lines) + "\n")
                if input_backs:
                    for input_back in input_backs:
                        with open(input_back) as in_handle:
                            for line in in_handle:
                                if len(line.split()) == 2:
                                    out_handle.write(line)
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

"""
Seq2c columns

Sample
Gene
Chromosome
Start
End
Length
Log2ratio Gene level log2 ratio (median)
Sig Significance (p-value if >= 3 exons, MAD value if < 3 exons)
BP_Whole BP or Whole: Value "Whole" means the gene is Del/Amp in whole, Value "BP" means with breakpoint. Empty means no copy number aberration
Amp_Del Amp or Del: Possible values are: "Amp", "Del", "NA", or empty
Ab_Seg: No. of segments/exons show aberration.
Total_Seg: Total segments/exons for the gene
Ab_log2ratio: The mean log2 ratio for exons with aberration.
Log2r_Diff Log2ratio Diff: The absolute log2 ratio difference between exons with aberration and other exons
Ab_Seg_Loc Exons and their location with aberration
Ab_Samples: No. of samples with the same exon(s) aberration (to catch false positives for exon with high variances)
Ab_Samples_Pcnt: Percentage of samples with the same exon(s) aberration
"""
