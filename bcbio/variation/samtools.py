"""Variant calling using samtools mpileup and bcftools.

http://samtools.sourceforge.net/mpileup.shtml
"""
import os
import subprocess

import sh

from bcbio import broad
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.variation.genotype import write_empty_vcf

def run_samtools(align_bam, ref_file, config, dbsnp=None, region=None,
                 out_file=None):
    """Detect SNPs and indels with samtools mpileup and bcftools.
    """
    broad_runner = broad.runner_from_config(config)
    broad_runner.run_fn("picard_index", align_bam)
    if out_file is None:
        out_file = "%s-variants.vcf" % os.path.splitext(align_bam)[0]
    if not file_exists(out_file):
        logger.info("Genotyping with samtools: {region} {fname}".format(
            region=region, fname=os.path.basename(align_bam)))
        variant_regions = config["algorithm"].get("variant_regions", None)
        target_regions = subset_variant_regions(variant_regions, region, out_file)
        if variant_regions is not None and not os.path.isfile(target_regions):
            write_empty_vcf(out_file)
        else:
            with file_transaction(out_file) as tx_out_file:
                _call_variants_samtools(align_bam, ref_file, config, target_regions,
                                        tx_out_file)
    return out_file

def _call_variants_samtools(align_bam, ref_file, config, target_regions, out_file):
    """Call variants with samtools in target_regions.
    """
    with open(out_file, "w") as out_handle:
        mpileup = sh.samtools.mpileup.bake(align_bam,
                                           f=ref_file, d=1000, L=1000,
                                           m=3, F=0.0002, D=True, S=True, u=True)
        if target_regions:
            mpileup=mpileup.bake(l=target_regions)
        bcftools = sh.bcftools.view.bake("-", v=True, c=True, g=True,
                                         _out=out_handle)
        bcftools(mpileup())
    
