"""Bayesian variant calling with FreeBayes.

http://bioinformatics.bc.edu/marthlab/FreeBayes
"""
import os
import shutil

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.provenance import do
from bcbio.variation import annotation, ploidy
from bcbio.variation.vcfutils import get_paired_bams, is_paired_analysis


def region_to_freebayes(region):
    if isinstance(region, (list, tuple)):
        chrom, start, end = region
        return "%s:%s..%s" % (chrom, start, end)
    else:
        return region

def _freebayes_options_from_config(items, config, out_file, region=None):
    opts = []
    opts += ["--ploidy", str(ploidy.get_ploidy(items, region))]

    variant_regions = utils.get_in(config, ("algorithm", "variant_regions"))
    target = subset_variant_regions(variant_regions, region, out_file)
    if target:
        if isinstance(target, basestring) and os.path.isfile(target):
            opts += ["--targets", target]
        else:
            opts += ["--region", region_to_freebayes(target)]
    resources = config_utils.get_resources("freebayes", config)
    if resources.get("options"):
        opts += resources["options"]
    return opts

def run_freebayes(align_bams, items, ref_file, assoc_files, region=None,
                  out_file=None):
    """Run FreeBayes variant calling, either paired tumor/normal or germline calling.
    """
    if is_paired_analysis(align_bams, items):
        call_file = _run_freebayes_paired(align_bams, items, ref_file,
                                          assoc_files, region, out_file)
    else:
        call_file = _run_freebayes_caller(align_bams, items, ref_file,
                                          assoc_files, region, out_file)

    return call_file

def _run_freebayes_caller(align_bams, items, ref_file, assoc_files,
                          region=None, out_file=None):
    """Detect SNPs and indels with FreeBayes.

    Performs post-filtering to remove very low quality variants which
    can cause issues feeding into GATK. Breaks variants into individual
    allelic primitives for analysis and evaluation.
    """
    config = items[0]["config"]
    if out_file is None:
        out_file = "%s-variants.vcf.gz" % os.path.splitext(align_bams[0])[0]
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            for align_bam in align_bams:
                bam.index(align_bam, config)
            freebayes = config_utils.get_program("freebayes", config)
            vcffilter = config_utils.get_program("vcffilter", config)
            vcfallelicprimitives = config_utils.get_program("vcfallelicprimitives", config)
            vcfstreamsort = config_utils.get_program("vcfstreamsort", config)
            input_bams = " ".join("-b %s" % x for x in align_bams)
            opts = " ".join(_freebayes_options_from_config(items, config, out_file, region))
            compress_cmd = "| bgzip -c" if out_file.endswith("gz") else ""
            cmd = ("{freebayes} -f {ref_file} {input_bams} {opts} | "
                   "{vcffilter} -f 'QUAL > 5' -s | {vcfallelicprimitives} | {vcfstreamsort} "
                   "{compress_cmd} > {tx_out_file}")
            do.run(cmd.format(**locals()), "Genotyping with FreeBayes", {})
    ann_file = annotation.annotate_nongatk_vcf(out_file, align_bams,
                                               assoc_files["dbsnp"],
                                               ref_file, config)
    return ann_file

def _run_freebayes_paired(align_bams, items, ref_file, assoc_files,
                          region=None, out_file=None):
    """Detect SNPs and indels with FreeBayes.

    This is used for paired tumor / normal samples.
    """
    config = items[0]["config"]
    if out_file is None:
        out_file = "%s-paired-variants.vcf.gz" % os.path.splitext(align_bams[0])[0]
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            paired = get_paired_bams(align_bams, items)
            if not paired.normal_bam:
                raise ValueError("Require both tumor and normal BAM files for FreeBayes cancer calling")

            vcfsamplediff = config_utils.get_program("vcfsamplediff", config)
            vcffilter = config_utils.get_program("vcffilter", config)
            freebayes = config_utils.get_program("freebayes", config)
            opts = " ".join(_freebayes_options_from_config(items, config, out_file, region))
            opts += " -f {}".format(ref_file)
            # NOTE: The first sample name in the vcfsamplediff call is
            # the one supposed to be the *germline* one

            # NOTE: -s in vcfsamplediff (strict checking: i.e., require no
            # reads in the germline to call somatic) is not used as it is
            # too stringent
            compress_cmd = "| bgzip -c" if out_file.endswith("gz") else ""
            cl = ("{freebayes} --pooled-discrete --genotype-qualities "
                  "{opts} {paired.tumor_bam} {paired.normal_bam} "
                  "| {vcffilter} -f 'QUAL > 1' -s "
                  "| {vcfsamplediff} VT {paired.normal_name} {paired.tumor_name} - "
                  "{compress_cmd} >  {tx_out_file}")
            bam.index(paired.tumor_bam, config)
            bam.index(paired.normal_bam, config)
            do.run(cl.format(**locals()), "Genotyping paired variants with FreeBayes", {})
    ann_file = annotation.annotate_nongatk_vcf(out_file, align_bams,
                                               assoc_files["dbsnp"], ref_file,
                                               config)
    return ann_file

def _move_vcf(orig_file, new_file):
    """Move a VCF file with associated index.
    """
    for ext in ["", ".idx"]:
        to_move = orig_file + ext
        if os.path.exists(to_move):
            shutil.move(to_move, new_file + ext)

def _clean_freebayes_output(line):
    """Clean FreeBayes output to make post-processing with GATK happy.

    XXX Not applied on recent versions which fix issues to be more compatible
    with bgzip output, but retained in case of need.

    - Remove lines from FreeBayes outputs where REF/ALT are identical:
      2       22816178        .       G       G       0.0339196
      or there are multiple duplicate alleles:
      4       60594753        .       TGAAA   T,T
    - Remove Type=Int specifications which are not valid VCF and GATK chokes
      on.
    """
    if line.startswith("#"):
        line = line.replace("Type=Int,D", "Type=Integer,D")
        return line
    else:
        parts = line.split("\t")
        alleles = [x.strip() for x in parts[4].split(",")] + [parts[3].strip()]
        if len(alleles) == len(set(alleles)):
            return line
    return None

def clean_vcf_output(orig_file, clean_fn, name="clean"):
    """Provide framework to clean a file in-place, with the specified clean
    function.
    """
    base, ext = utils.splitext_plus(orig_file)
    out_file = "{0}-{1}{2}".format(base, name, ext)
    if not utils.file_exists(out_file):
        with open(orig_file) as in_handle:
            with file_transaction(out_file) as tx_out_file:
                with open(tx_out_file, "w") as out_handle:
                    for line in in_handle:
                        update_line = clean_fn(line)
                        if update_line:
                            out_handle.write(update_line)
        _move_vcf(orig_file, "{0}.orig".format(orig_file))
        _move_vcf(out_file, orig_file)
        with open(out_file, "w") as out_handle:
            out_handle.write("Moved to {0}".format(orig_file))
