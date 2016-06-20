"""Cohort based copy number calling in gene regions using Seq2C.

Seq2C calls across multiple samples without explicit background samples,
using gene regions as segments.

This requires coverage calculation in each sample and gene, followed by global
calling across all samples.
"""
import os
import subprocess

from bcbio import utils
from bcbio.pipeline import datadict as dd
from bcbio.distributed.multi import run_multicore
from bcbio.distributed.transaction import file_transaction
from bcbio.structural import cnvkit

def precall(items):
    """Perform initial pre-calling steps -- coverage calcuation by sample.

    Coverage calculation per gene, using existing CNVkit coverage calculation
    parallelization in place of seq2cov.pl.
    """
    items = [utils.to_single_data(x) for x in items]
    assert len(items) == 1, "Expect one item to Seq2C coverage calculation"
    data = items[0]
    assert dd.get_coverage_interval(data) != "genome", "Seq2C only for amplicon and exome sequencing"
    work_dir = _sv_workdir(data)
    parallel = {"type": "local", "cores": dd.get_cores(data), "progs": ["seq2c"]}
    split_cnns = run_multicore(cnvkit._cnvkit_coverage,
                                [(dd.get_align_bam(data), bed, "", work_dir, data)
                                 for bed in cnvkit._split_bed(dd.get_variant_regions(data), data)],
                               data["config"], parallel)
    coverage_cnns = cnvkit._merge_coverage(split_cnns, data)
    assert len(coverage_cnns) == 1
    if "sv" not in data:
        data["sv"] = []
    data["sv"].append({"variantcaller": "seq2c",
                       "cnn": coverage_cnns[0]["file"]})
    return [data]

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
