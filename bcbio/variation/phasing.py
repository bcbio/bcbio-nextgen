"""Approaches for calculating haplotype phasing of variants.

Currently supports GATK's Read-Backed phasing:

http://www.broadinstitute.org/gsa/wiki/index.php/Read-backed_phasing_algorithm
"""
import os

from bcbio import broad
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction

def read_backed_phasing(vcf_file, bam_file, genome_file, config):
    """Annotate predicted variant effects using snpEff.
    """
    broad_runner = broad.runner_from_config(config)
    out_file = "%s-phased%s" % os.path.splitext(vcf_file)
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "ReadBackedPhasing",
                      "-R", genome_file,
                      "-I", bam_file,
                      "--variant", vcf_file,
                      "--out", tx_out_file]
            broad_runner.run_gatk(params)
    return out_file
