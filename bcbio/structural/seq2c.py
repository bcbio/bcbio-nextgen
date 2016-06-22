"""Cohort based copy number calling in gene regions using Seq2C.

Seq2C calls across multiple samples without explicit background samples,
using gene regions as segments.

This requires coverage calculation in each sample and gene, followed by global
calling across all samples.
"""
import os
import subprocess
from collections import defaultdict

from bcbio import utils
from bcbio.pipeline import datadict as dd
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.log import logger
from pipeline import config_utils


def precall(items):
    """Perform initial pre-calling steps -- coverage calcuation by sample.

    Use sambamba to call average region coverage in regions
    """
    items = [utils.to_single_data(x) for x in items]
    assert len(items) == 1, "Expect one item to Seq2C coverage calculation"
    data = items[0]
    assert dd.get_coverage_interval(data) != "genome", "Seq2C only for amplicon and exome sequencing"

    bed_file = dd.get_variant_regions(data)
    bam_file = dd.get_align_bam(data)
    sample_name = data['name'][1]

    work_dir = _sv_workdir(data)
    sambamba_depth_file = os.path.join(work_dir, sample_name + '.sambamba_depth.txt')
    sambamba = config_utils.get_program("sambamba", data["config"])
    num_cores = dd.get_cores(data)
    if not utils.file_exists(sambamba_depth_file):
        with file_transaction(data, sambamba_depth_file) as tx_out_file:
            cmd = "{sambamba} depth region -t {num_cores} -L {bed_file} -o {tx_out_file} {bam_file}"
            do.run(cmd.format(**locals()), "Calling sambamba region depth")
    logger.debug('Saved to ' + sambamba_depth_file)

    logger.debug('Converting sambamba depth output to cov2lr.pl input')
    out_file = os.path.join(work_dir, os.path.splitext(os.path.basename(bam_file))[0] + '.seq2cov.txt')
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            sambabma_depth_to_seq2cov(sambamba_depth_file, tx_out_file, sample_name)
    logger.debug('Saved to ' + out_file)

    if "sv" not in data:
        data["sv"] = []
    data["sv"].append({"variantcaller": "seq2c",
                       "cov": out_file})
    return [data]

''' sambamba_depth_output_fpath:
# chrom chromStart  chromEnd  F3       readCount  minDepth  meanCoverage  stdDev   percentage1  percentage5  percentage10  ...  sampleName
chr20   68345       68413     DEFB125  56         28        32.5          1.66716  100          100          100           ...  chr20_tumor
chr20   76640       77301     DEFB125  279        24        36.9213       5.74231  100          100          100           ...  chr20_tumor
'''
''' seq2cov:
chr20_tumor_1   DEFB125   chr20   68346   68413   Amplicon    68   28.0
chr20_tumor_1   DEFB125   chr20   76641   77301   Amplicon    661  24.0
chr20_tumor_1   DEFB125   chr20   68346   77301   Whole-Gene  729  24.3731138546
chr20_tumor_1   DEFB126   chr20   123247  123332  Amplicon    86   40.0
'''
def sambabma_depth_to_seq2cov(input_fpath, output_fpath, sample_name):
    with open(input_fpath) as f:
        ave_depth_col = next(f).split('\t').index('meanCoverage')

    # first round: collecting gene ends
    gene_end_by_gene = defaultdict(lambda: -1)
    with open(input_fpath) as f:
        for l in f:
            if l.startswith('#'): continue
            fs = l.replace('\n', '').split('\t')
            if any(fs[i] == '.' for i in [0, 1, 2, 3, ave_depth_col]): continue
            end = int(fs[2])
            gene_name = fs[3]
            gene_end_by_gene[gene_name] = max(gene_end_by_gene[gene_name], end)

    # second round: calculating coverage
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

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "seq2c"))

def run(items, background=None):
    """Normalization and log2 ratio calculation plus CNV calling for full cohort.

    - Prepare coverage file in correct format
    - Prepare read counts for each sample
    - cov2lr.pl -- log2 ratio calculation (do we need this with CNNs from CNVkit?)
    - lr2gene.pl -- call amplifications and deletions
    """
    items = [utils.to_single_data(x) for x in items]
    work_dir = _sv_workdir(items[0])
    coverage_file = _combine_coverages(items, work_dir)
    read_mapping_file = _calculate_mapping_reads(items, work_dir)
    return items

def _combine_coverages(items, work_dir):
    """Combine coverage cnns calculated for individual inputs into single file.
    """
    out_file = os.path.join(work_dir, "sample_coverages.txt")
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            for data in items:
                svouts = [x for x in data["sv"] if x["variantcaller"] == "seq2c"]
                assert len(svouts) == 1
                cnn_file = svouts[0]["cnn"]
                print cnn_file
    return out_file

def _calculate_mapping_reads(items, work_dir):
    """Calculate read counts from samtools idxstats for each sample.
    """
    out_file = os.path.join(work_dir, "mapping_reads.txt")
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                for data in items:
                    count = 0
                    for line in subprocess.check_output(["samtools", "idxstats",
                                                         dd.get_align_bam(data)]).split("\n"):
                        if line.strip():
                            count += int(line.split("\t")[2])
                    out_handle.write("%s\t%s\n" % (dd.get_sample_name(data), count))
    return out_file
