"""Variant calling using samtools mpileup and bcftools.

http://samtools.sourceforge.net/mpileup.shtml
"""
import os

import sh

from bcbio import broad
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.variation.genotype import write_empty_vcf
from bcbio.variation.realign import has_aligned_reads

def shared_variantcall(call_fn, name, align_bams, ref_file, config,
                       dbsnp=None, region=None, out_file=None):
    """Provide base functionality for prepping and indexing for variant calling.
    """
    broad_runner = broad.runner_from_config(config)
    for x in align_bams:
        broad_runner.run_fn("picard_index", x)
    if out_file is None:
        out_file = "%s-variants.vcf" % os.path.splitext(align_bams[0])[0]
    if not file_exists(out_file):
        logger.info("Genotyping with {name}: {region} {fname}".format(name=name,
            region=region, fname=os.path.basename(align_bams[0])))
        variant_regions = config["algorithm"].get("variant_regions", None)
        target_regions = subset_variant_regions(variant_regions, region, out_file)
        if ((variant_regions is not None and not os.path.isfile(target_regions))
              or not all(has_aligned_reads(x, region) for x in align_bams)):
            write_empty_vcf(out_file)
        else:
            with file_transaction(out_file) as tx_out_file:
                call_fn(align_bams, ref_file, config, target_regions,
                        tx_out_file)
    return out_file


def run_samtools(align_bams, ref_file, config, dbsnp=None, region=None,
                 out_file=None):
    """Detect SNPs and indels with samtools mpileup and bcftools.
    """
    return shared_variantcall(_call_variants_samtools, "samtools", align_bams, ref_file,
                              config, dbsnp, region, out_file)

def prep_mpileup(align_bams, ref_file, max_read_depth, target_regions=None, want_bcf=True):
    mpileup = sh.samtools.mpileup.bake(*align_bams,
                                       f=ref_file, d=max_read_depth, L=max_read_depth,
                                       m=3, F=0.0002)
    if want_bcf:
        mpileup = mpileup.bake(D=True, S=True, u=True)
    if target_regions:
        mpileup = mpileup.bake(l=target_regions)
    return mpileup

def _call_variants_samtools(align_bams, ref_file, config, target_regions, out_file):
    """Call variants with samtools in target_regions.
    """
    max_read_depth = 1000
    with open(out_file, "w") as out_handle:
        mpileup = prep_mpileup(align_bams, ref_file, max_read_depth, target_regions)
        bcftools = sh.bcftools.view.bake("-", v=True, c=True, g=True)
        varfilter = sh.Command("vcfutils.pl").varFilter.bake(D=max_read_depth, _out=out_handle)
        varfilter(bcftools(mpileup()))
