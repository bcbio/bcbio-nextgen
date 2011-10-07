"""Annotated variant VCF files with additional information.

- GATK variant annotation with snpEff predicted effects.
"""
import os

from bcbio import broad
from bcbio.utils import file_transaction

# ## snpEff annotation

def annotate_effects(orig_file, snpeff_file, genome_file, config):
    """Annotate predicted variant effects using snpEff.
    """
    broad_runner = broad.runner_from_config(config)
    out_file = "%s-annotated%s" % os.path.splitext(orig_file)
    params = ["-T", "VariantAnnotator",
              "-R", genome_file,
              "-A", "SnpEff",
              "--variant", orig_file,
              "--snpEffFile", snpeff_file,
              "--out", out_file]
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        with file_transaction(out_file):
            broad_runner.run_gatk(params)
    return out_file
