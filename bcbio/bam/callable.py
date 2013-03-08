"""Examine callable regions following genome mapping of short reads.

Identifies callable analysis regions surrounded by larger regions lacking
aligned bases. This allows parallelization of smaller chromosome chunks
through post-processing and variant calling, with each sub-section
mapping handled separately.
"""
import os

from bcbio import utils, broad
from bcbio.log import logger
from bcbio.distributed.transaction import file_transaction

def calc_callable_loci(in_bam, ref_file, config):
    """Determine callable bases for input BAM using Broad's CallableLoci walker.
    """
    broad_runner = broad.runner_from_config(config)
    out_bed = "%s-callable.bed" % os.path.splitext(in_bam)[0]
    out_summary = "%s-callable-summary.txt" % os.path.splitext(in_bam)[0]
    variant_regions = config["algorithm"].get("variant_regions", None)
    if not utils.file_exists(out_bed):
        with file_transaction(out_bed) as tx_out_file:
            broad_runner.run_fn("picard_index", in_bam)
            params = ["-T", "CallableLoci",
                      "-R", ref_file,
                      "-I", in_bam,
                      "--out", tx_out_file,
                      "--summary", out_summary]
            if variant_regions:
                params += ["-L", variant_regions]
            broad_runner.run_gatk(params)
    return out_bed

def block_regions(in_bam, ref_file, config):
    """Find blocks of regions for analysis from mapped input BAM file.
    """
    callable_bed = calc_callable_loci(in_bam, ref_file, config)
    return []
