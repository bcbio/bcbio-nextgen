"""Variant calling using samtools mpileup and bcftools.

http://samtools.sourceforge.net/mpileup.shtml
"""
import os

from bcbio import broad
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.provenance import do
from bcbio.variation import annotation, bamprep, realign, vcfutils

def shared_variantcall(call_fn, name, align_bams, ref_file, config,
                       assoc_files, region=None, out_file=None):
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
        if ((variant_regions is not None and isinstance(target_regions, basestring)
              and not os.path.isfile(target_regions))
              or not all(realign.has_aligned_reads(x, region) for x in align_bams)):
            vcfutils.write_empty_vcf(out_file)
        else:
            with file_transaction(out_file) as tx_out_file:
                call_fn(align_bams, ref_file, config, target_regions,
                        tx_out_file)
    ann_file = annotation.annotate_nongatk_vcf(out_file, align_bams, assoc_files.dbsnp,
                                               ref_file, config)
    return ann_file


def run_samtools(align_bams, items, ref_file, assoc_files, region=None,
                 out_file=None):
    """Detect SNPs and indels with samtools mpileup and bcftools.
    """
    return shared_variantcall(_call_variants_samtools, "samtools", align_bams, ref_file,
                              items[0]["config"], assoc_files, region, out_file)

def prep_mpileup(align_bams, ref_file, max_read_depth, config,
                 target_regions=None, want_bcf=True):
    cl = [config_utils.get_program("samtools", config), "mpileup",
          "-f", ref_file, "-d", str(max_read_depth), "-L", str(max_read_depth),
          "-m", "3", "-F", "0.0002"]
    if want_bcf:
        cl += ["-D", "-S", "-u"]
    if target_regions:
        str_regions = bamprep.region_to_gatk(target_regions)
        if os.path.isfile(str_regions):
            cl += ["-l", str_regions]
        else:
            cl += ["-r", str_regions]
    cl += align_bams
    return " ".join(cl)

def _call_variants_samtools(align_bams, ref_file, config, target_regions, out_file):
    """Call variants with samtools in target_regions.
    """
    max_read_depth = "1000"
    mpileup = prep_mpileup(align_bams, ref_file, max_read_depth, config,
                           target_regions=target_regions)
    bcftools = config_utils.get_program("bcftools", config)
    vcfutils = config_utils.get_program("vcfutils.pl", config)
    cmd = ("{mpileup} "
           "| {bcftools} view -v -c -g - "
           "| {vcfutils} varFilter -D {max_read_depth} "
           "> {out_file}")
    logger.info(cmd.format(**locals()))
    do.run(cmd.format(**locals()), "Variant calling with samtools", {})
