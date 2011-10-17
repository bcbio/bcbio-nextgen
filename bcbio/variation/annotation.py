"""Annotated variant VCF files with additional information.

- GATK variant annotation with snpEff predicted effects.
"""
import os

from bcbio import broad
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction

# ## snpEff annotation

def annotate_effects(orig_file, snpeff_file, genome_file, config):
    """Annotate predicted variant effects using snpEff.
    """
    broad_runner = broad.runner_from_config(config)
    out_file = "%s-annotated%s" % os.path.splitext(orig_file)
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "VariantAnnotator",
                      "-R", genome_file,
                      "-A", "SnpEff",
                      "--variant", orig_file,
                      "--snpEffFile", snpeff_file,
                      "--out", tx_out_file]
            broad_runner.run_gatk(params)
    return out_file
