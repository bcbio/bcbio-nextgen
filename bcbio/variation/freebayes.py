"""Bayesian variant calling with FreeBayes.

http://bioinformatics.bc.edu/marthlab/FreeBayes
"""
import os
import shutil

from bcbio import broad
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.provenance import do
from bcbio.variation import annotation

def region_to_freebayes(region):
    if isinstance(region, (list, tuple)):
        chrom, start, end = region
        return "%s:%s..%s" % (chrom, start, end)
    else:
        return region

def _freebayes_options_from_config(aconfig, out_file, region=None):
    opts = []
    ploidy = aconfig.get("ploidy", 2)
    opts += ["--ploidy", str(ploidy)]

    variant_regions = aconfig.get("variant_regions", None)
    target = subset_variant_regions(variant_regions, region, out_file)
    if target:
        opts += ["--region" if target == region else "--targets",
                 region_to_freebayes(target)]
    #background = aconfig.get("call_background", None)
    #if background and os.path.exists(background):
    #    opts += ["--variant-input", background]
    return opts

def run_freebayes(align_bams, items, ref_file, assoc_files, region=None,
                  out_file=None):
    """Detect SNPs and indels with FreeBayes.
    """
    config = items[0]["config"]
    broad_runner = broad.runner_from_config(config)
    if out_file is None:
        out_file = "%s-variants.vcf" % os.path.splitext(align_bams[0])[0]
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            cl = [config_utils.get_program("freebayes", config),
                  "-v", tx_out_file, "-f", ref_file,
                  "--use-mapping-quality", "--pvar", "0.7"]
            for align_bam in align_bams:
                broad_runner.run_fn("picard_index", align_bam)
                cl += ["-b", align_bam]
            cl += _freebayes_options_from_config(config["algorithm"], out_file, region)
            do.run(cl, "Genotyping with FreeBayes", {})
        _clean_freebayes_output(out_file)
    ann_file = annotation.annotate_nongatk_vcf(out_file, align_bams, assoc_files.dbsnp,
                                               ref_file, config)
    return ann_file

def _move_vcf(orig_file, new_file):
    """Move a VCF file with associated index.
    """
    for ext in ["", ".idx"]:
        to_move = orig_file + ext
        if os.path.exists(to_move):
            shutil.move(to_move, new_file + ext)

def _clean_freebayes_output(in_file):
    """Clean FreeBayes output to make post-processing with GATK happy.
    - Remove lines from FreeBayes outputs where REF/ALT are identical:
      2       22816178        .       G       G       0.0339196
    - Remove Type=Int specifications which are not valid VCF and GATK chokes on.
    """
    out_file = apply("{0}-nodups{1}".format, os.path.splitext(in_file))
    if not file_exists(out_file):
        with open(in_file) as in_handle:
            with open(out_file, "w") as out_handle:
                for line in in_handle:
                    if line.startswith("#"):
                        line = line.replace("Type=Int,D", "Type=Integer,D")
                        out_handle.write(line)
                    else:
                        parts = line.split("\t")
                        if parts[3] != parts[4]:
                            out_handle.write(line)
        _move_vcf(in_file, "{0}.orig".format(in_file))
        _move_vcf(out_file, in_file)
        with open(out_file, "w") as out_handle:
            out_handle.write("Moved to {0}".format(in_file))
