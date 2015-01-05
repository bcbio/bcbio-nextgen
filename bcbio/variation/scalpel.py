"""InDel calling using Scalpel

https://sourceforge.net/p/scalpel/code/ci/master/tree/
"""

from __future__ import print_function
import os

try:
    import vcf
except ImportError:
    vcf = None

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import subset_variant_regions, remove_lcr_regions
from bcbio.provenance import do
from bcbio.variation import annotation
from bcbio.variation.vcfutils import get_paired_bams, is_paired_analysis, bgzip_and_index

def _scalpel_options_from_config(items, config, out_file, region, tmp_path):
    opts = []
    # output vcf, report only variants within bed regions
    opts += ["--format", "vcf", "--intarget"]
    # Improve sensitivity in low coverage regions
    opts += ["--covthr 3", "--lowcov 1"]
    # Avoid oversampling in repeat regions
    opts += ["--pathlimit", "10000"]
    variant_regions = utils.get_in(config, ("algorithm", "variant_regions"))
    target = subset_variant_regions(variant_regions, region, out_file, items)
    if target:
        if isinstance(target, basestring) and os.path.isfile(target):
            target_bed = target
        else:
            target_bed = os.path.join(tmp_path, "tmp.bed")
            with file_transaction(config, target_bed) as tx_tmp_bed:
                if not isinstance(region, (list, tuple)):
                    message = ("Region must be a tuple - something odd just happened")
                    raise ValueError(message)
                chrom, start, end = region
                with open(tx_tmp_bed, "w") as out_handle:
                    print("%s\t%s\t%s" % (chrom, start, end), file=out_handle)
        opts += ["--bed", remove_lcr_regions(target_bed, items)]
    resources = config_utils.get_resources("scalpel", config)
    if resources.get("options"):
        opts += resources["options"]
    if "--outratio" not in " ".join(opts):
        # add minimum reportable allele frequency, for which Scalpel defaults to 5
        # but other somatic tools in bcbio default to 10
        min_af = float(utils.get_in(config, ("algorithm",
                                             "min_allele_fraction"), 10)) / 100.0
        opts += ["--outratio", str(min_af)]
    return opts

def is_installed(config):
    """Check for scalpel installation on machine.
    """
    try:
        config_utils.get_program("scalpel", config)
        return True
    except config_utils.CmdNotFound:
        return False

def run_scalpel(align_bams, items, ref_file, assoc_files, region=None,
                  out_file=None):
    """Run Scalpel indel calling, either paired tumor/normal or germline calling.
    """
    if region is None:
        message = ("A region must be provided for Scalpel")
        raise ValueError(message)
    if is_paired_analysis(align_bams, items):
        call_file = _run_scalpel_paired(align_bams, items, ref_file,
                                          assoc_files, region, out_file)
    else:
        call_file = _run_scalpel_caller(align_bams, items, ref_file,
                                          assoc_files, region, out_file)
    return call_file

def _run_scalpel_caller(align_bams, items, ref_file, assoc_files,
                          region=None, out_file=None):
    """Detect indels with Scalpel.

    Single sample mode.
    """
    config = items[0]["config"]
    if out_file is None:
        out_file = "%s-variants.vcf.gz" % os.path.splitext(align_bams[0])[0]
    if not utils.file_exists(out_file):
        with file_transaction(config, out_file) as tx_out_file:
            for align_bam in align_bams:
                bam.index(align_bam, config)
            scalpel = config_utils.get_program("scalpel", config)
            if len(align_bams) > 1:
                message = ("Scalpel does not currently support batch calling!")
                raise ValueError(message)
            input_bams = " ".join("%s" % x for x in align_bams)
            tmp_path = os.path.dirname(tx_out_file)
            opts = " ".join(_scalpel_options_from_config(items, config, out_file, region, tmp_path))
            opts += " --dir %s" % tmp_path
            min_cov = "3"  # minimum coverage
            opts += " --mincov %s" % min_cov
            cmd = ("{scalpel} --single {opts} --ref {ref_file} --bam {input_bams} ")
            # first run into temp folder
            do.run(cmd.format(**locals()), "Genotyping with Scalpel", {})
            # parse produced variant file further
            scalpel_tmp_file = bgzip_and_index(os.path.join(tmp_path, "variants." + min_cov + "x.indel.vcf"), config)
            compress_cmd = "| bgzip -c" if out_file.endswith("gz") else ""
            bcftools_cmd_chi2 = get_scalpel_bcftools_filter_expression("chi2", config)
            sample_name_str = items[0]["name"][1]
            cl2 = ("{bcftools_cmd_chi2} {scalpel_tmp_file} | sed 's/sample_name/{sample_name_str}/g' | "
                   "vcfallelicprimitives --keep-geno | vcffixup - | vcfstreamsort "
                   "{compress_cmd} > {tx_out_file}")
            do.run(cl2.format(**locals()), "Finalising Scalpel variants", {})
    ann_file = annotation.annotate_nongatk_vcf(out_file, align_bams,
                                               assoc_files.get("dbsnp"),
                                               ref_file, config)
    return ann_file

def _run_scalpel_paired(align_bams, items, ref_file, assoc_files,
                          region=None, out_file=None):
    """Detect indels with Scalpel.

    This is used for paired tumor / normal samples.
    """
    config = items[0]["config"]
    if out_file is None:
        out_file = "%s-paired-variants.vcf.gz" % os.path.splitext(align_bams[0])[0]
    if not utils.file_exists(out_file):
        with file_transaction(config, out_file) as tx_out_file:
            paired = get_paired_bams(align_bams, items)
            if not paired.normal_bam:
                ann_file = _run_scalpel_caller(align_bams, items, ref_file,
                                               assoc_files, region, out_file)
                return ann_file
            vcffilter = config_utils.get_program("vcffilter", config)
            scalpel = config_utils.get_program("scalpel", config)
            vcfstreamsort = config_utils.get_program("vcfstreamsort", config)
            tmp_path = os.path.dirname(tx_out_file)
            opts = " ".join(_scalpel_options_from_config(items, config, out_file, region, tmp_path))
            opts += " --ref {}".format(ref_file)
            opts += " --dir %s" % tmp_path
            min_cov = "3"  # minimum coverage
            opts += " --mincov %s" % min_cov
            cl = ("{scalpel} --somatic {opts} --tumor {paired.tumor_bam} --normal {paired.normal_bam}")
            bam.index(paired.tumor_bam, config)
            bam.index(paired.normal_bam, config)
            do.run(cl.format(**locals()), "Genotyping paired variants with Scalpel", {})
            # somatic
            scalpel_tmp_file = bgzip_and_index(os.path.join(tmp_path, "main/somatic." + min_cov + "x.indel.vcf"), config)
            # common
            scalpel_tmp_file_common = bgzip_and_index(os.path.join(tmp_path, "main/common." + min_cov + "x.indel.vcf"), config)
            compress_cmd = "| bgzip -c" if out_file.endswith("gz") else ""
            bcftools_cmd_chi2 = get_scalpel_bcftools_filter_expression("chi2", config)
            bcftools_cmd_common = get_scalpel_bcftools_filter_expression("reject", config)
            cl2 = ("vcfcat <({bcftools_cmd_chi2} {scalpel_tmp_file}) "
                   "<({bcftools_cmd_common} {scalpel_tmp_file_common}) | "
                   " sed 's/sample_name/{paired.tumor_name}/g' | "
                   "{vcfstreamsort} {compress_cmd} > {tx_out_file}")
            do.run(cl2.format(**locals()), "Finalising Scalpel variants", {})
            
    ann_file = annotation.annotate_nongatk_vcf(out_file, align_bams,
                                               assoc_files.get("dbsnp"), ref_file,
                                               config)
    return ann_file

def get_scalpel_bcftools_filter_expression(filter_type, config):
    bcftools = config_utils.get_program("bcftools", config)
    filter_string = "{bcftools} filter -m '+' -O v --soft-filter "
    if filter_type == "chi2":
        filter_string += "'CHI2FILTER' -e 'INFO/CHI2 > 20.0' "
    elif filter_type == "reject":
        filter_string += "'REJECT' -e '%TYPE=\"indel\"' "
    else:
        return "zcat"
    return filter_string.format(**locals())
