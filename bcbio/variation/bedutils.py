"""Utilities for manipulating BED files.
"""
import os
import shutil

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.variation import vcfutils

def clean_file(in_file, data, prefix=""):
    """Prepare a clean input BED file without headers or overlapping segments.

    Overlapping regions (1:1-100, 1:90-100) cause issues with callers like FreeBayes
    that don't collapse BEDs prior to using them.
    """
    bedtools = config_utils.get_program("bedtools", data["config"])
    if in_file:
        bedprep_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "bedprep"))
        out_file = os.path.join(bedprep_dir, "%s%s" % (prefix, os.path.basename(in_file)))
        if not utils.file_exists(out_file):
            with file_transaction(out_file) as tx_out_file:
                cmd = "sort -k1,1 -k2,2n {in_file} | {bedtools} merge -i > {tx_out_file}"
                do.run(cmd.format(**locals()), "Prepare cleaned BED file", data)
        vcfutils.bgzip_and_index(out_file, data["config"], remove_orig=False)
        return out_file

def clean_inputs(data):
    """Clean BED input files to avoid overlapping segments that cause downstream issues.
    """
    data["config"]["algorithm"]["variant_regions"] = clean_file(
        utils.get_in(data, ("config", "algorithm", "variant_regions")), data)
    return data

def combine(in_files, out_file, config):
    """Combine multiple BED files into a single output.
    """
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                for in_file in in_files:
                    with open(in_file) as in_handle:
                        shutil.copyfileobj(in_handle, out_handle)
    return out_file
