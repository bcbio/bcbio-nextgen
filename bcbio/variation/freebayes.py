"""Bayesian variant calling with FreeBayes.

http://bioinformatics.bc.edu/marthlab/FreeBayes
"""
import os
import subprocess

from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import log

def run_freebayes(align_bam, ref_file, config, dbsnp=None, region=None,
                  out_file=None):
    """Detect small polymorphisms with FreeBayes.
    """
    if out_file is None:
        out_file = "%s-variants.vcf" % os.path.splitext(align_bam)[0]
    if not file_exists(out_file):
        log.info("Genotyping with FreeBayes: {region} {fname}".format(
            region=region,
            fname=os.path.basename(align_bam)))
        with file_transaction(out_file) as tx_out_file:
            cl = [config["program"].get("freebayes", "freebayes"),
                  "-b", align_bam, "-v", tx_out_file, "-f", ref_file]
            if region:
                cl.extend(["-r", region])
            try:
                subprocess.check_call(cl)
            # XXX Temporary, work around freebayes issue; need to recall these regions
            # later so this is an ugly silent fix. Will need to grep for 'freebayes failed'
            # https://github.com/ekg/freebayes/issues/22
            except subprocess.CalledProcessError:
                with open(tx_out_file, "w") as out_handle:
                    out_handle.write("##fileformat=VCFv4.1\n"
                                     "## No variants; freebayes failed\n"
                                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    return out_file
