"""Variant calling using samtools 1.0 mpileup and bcftools.

http://www.htslib.org/workflow/#mapping_to_variant
"""
import os
from distutils.version import LooseVersion

import toolz as tz

from bcbio import bam
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.provenance import do, programs
from bcbio.variation import annotation, bamprep, bedutils, vcfutils

def shared_variantcall(call_fn, name, align_bams, ref_file, items,
                       assoc_files, region=None, out_file=None):
    """Provide base functionality for prepping and indexing for variant calling.
    """
    config = items[0]["config"]
    if out_file is None:
        if vcfutils.is_paired_analysis(align_bams, items):
            out_file = "%s-paired-variants.vcf.gz" % config["metdata"]["batch"]
        else:
            out_file = "%s-variants.vcf.gz" % os.path.splitext(align_bams[0])[0]
    if not file_exists(out_file):
        logger.debug("Genotyping with {name}: {region} {fname}".format(
              name=name, region=region, fname=os.path.basename(align_bams[0])))
        for x in align_bams:
            bam.index(x, config)
        variant_regions = bedutils.merge_overlaps(tz.get_in(["config", "algorithm", "variant_regions"], items[0]),
                                                  items[0])
        target_regions = subset_variant_regions(variant_regions, region, out_file)
        if (variant_regions is not None and isinstance(target_regions, basestring)
              and not os.path.isfile(target_regions)):
            vcfutils.write_empty_vcf(out_file, config)
        else:
            with file_transaction(config, out_file) as tx_out_file:
                call_fn(align_bams, ref_file, items, target_regions,
                        tx_out_file)
    if out_file.endswith(".gz"):
        out_file = vcfutils.bgzip_and_index(out_file, config)
    ann_file = annotation.annotate_nongatk_vcf(out_file, align_bams, assoc_files.get("dbsnp"),
                                               ref_file, config)
    return ann_file

def run_samtools(align_bams, items, ref_file, assoc_files, region=None,
                 out_file=None):
    """Detect SNPs and indels with samtools mpileup and bcftools.
    """
    return shared_variantcall(_call_variants_samtools, "samtools", align_bams, ref_file,
                              items, assoc_files, region, out_file)

def prep_mpileup(align_bams, ref_file, config, max_read_depth=None,
                 target_regions=None, want_bcf=True):
    cl = [config_utils.get_program("samtools", config), "mpileup", "-f", ref_file]
    if max_read_depth:
        cl += ["-d", str(max_read_depth), "-L", str(max_read_depth)]
    if want_bcf:
        cl += ["-t", "DP", "-u", "-g"]
    if target_regions:
        str_regions = bamprep.region_to_gatk(target_regions)
        if os.path.isfile(str_regions):
            cl += ["-l", str_regions]
        else:
            cl += ["-r", str_regions]
    cl += align_bams
    return " ".join(cl)

def _call_variants_samtools(align_bams, ref_file, items, target_regions, tx_out_file):
    """Call variants with samtools in target_regions.

    Works around a GATK VCF 4.2 compatibility issue in samtools 1.0
    by removing addition 4.2-only isms from VCF header lines.
    """
    config = items[0]["config"]
    mpileup = prep_mpileup(align_bams, ref_file, config,
                           target_regions=target_regions, want_bcf=True)
    bcftools = config_utils.get_program("bcftools", config)
    bcftools_version = programs.get_version("bcftools", config=config)
    samtools_version = programs.get_version("samtools", config=config)
    if LooseVersion(samtools_version) <= LooseVersion("0.1.19"):
        raise ValueError("samtools calling not supported with pre-1.0 samtools")
    bcftools_opts = "call -v -m"
    compress_cmd = "| bgzip -c" if tx_out_file.endswith(".gz") else ""
    cmd = ("{mpileup} "
           "| {bcftools} {bcftools_opts} - "
           "| vt normalize -n -q -r {ref_file} - "
           "| sed 's/VCFv4.2/VCFv4.1/' "
           "| sed 's/,Version=3>/>/' "
           "| sed 's/,Version=\"3\">/>/' "
           "| sed 's/Number=R/Number=./' "
           "{compress_cmd} > {tx_out_file}")
    do.run(cmd.format(**locals()), "Variant calling with samtools", items[0])
