"""InDel calling using Scalpel

https://sourceforge.net/p/scalpel/code/ci/master/tree/
"""

from __future__ import print_function
from collections import namedtuple
import os
import shutil

try:
    import vcf
except ImportError:
    vcf = None

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.provenance import do
from bcbio.variation import annotation, ploidy
from bcbio.variation.vcfutils import get_paired_bams, is_paired_analysis, bgzip_and_index

def _scalpel_options_from_config(items, config, out_file, region, tmp_path):
    opts = []
    opts += ["--format","vcf","--intarget"] # output vcf, report only variants within bed regions
    variant_regions = utils.get_in(config, ("algorithm", "variant_regions"))
    target = subset_variant_regions(variant_regions, region, out_file, items)
    if target:
        if isinstance(target, basestring) and os.path.isfile(target):
            opts += ["--bed", target]
        else:
            tmp_bed = os.path.join(tmp_path,"tmp.bed")
            with file_transaction(tmp_bed) as tx_tmp_bed:
                if not isinstance(region, (list, tuple)):
                    message = ("Region must be a tuple - something odd just happened")
                    raise ValueError(message)
                chrom, start, end = region
                print("%s\t%s\t%s" % (chrom, start, end), file=tx_tmp_bed)
            opts += ["--bed", tmp_bed]
    resources = config_utils.get_resources("scalpel", config)
    if resources.get("options"):
        opts += resources["options"]
    if "--outratio" not in " ".join(opts):
        # add minimum reportable allele frequency, for which Scalpel defaults to 5 but other somatic tools in bcbio default to 10
         min_af = float(utils.get_in(config, ("algorithm",
                                              "min_allele_fraction"),10)) / 100.0
         opts += ["--outratio", str(min_af)]
    return opts

def run_scalpel(align_bams, items, ref_file, assoc_files, region=None,
                  out_file=None):
    """Run Scalpel indel calling, either paired tumor/normal or germline calling.
    """
    if region == None:
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
        with file_transaction(out_file) as tx_out_file:
            for align_bam in align_bams:
                bam.index(align_bam, config)
            scalpel = config_utils.get_program("scalpel", config)
            vcfallelicprimitives = config_utils.get_program("vcfallelicprimitives", config)
            vcfstreamsort = config_utils.get_program("vcfstreamsort", config)
            if len(align_bams) > 1:
                message = ("Scalpel does not currently support batch calling!")
                raise ValueError(message)
            input_bams = " ".join("%s" % x for x in align_bams)
            tmp_path = os.path.dirname(tx_out_file)
            opts = " ".join(_scalpel_options_from_config(items, config, out_file, region, tmp_path))
            opts += " --dir %s" % tmp_path
            min_cov = "5" # minimum coverage (default 5)
            opts += " --mincov %s" % min_cov
            cmd = ("{scalpel} --single {opts} --ref {ref_file} --bam {input_bams} ")
            # first run into temp folder
            do.run(cmd.format(**locals()), "Genotyping with Scalpel", {})
            # parse produced variant file further
            scalpel_tmp_file = os.path.join(tmp_path, "variants."+min_cov+"x.indel.vcf")
            compress_cmd = "| bgzip -c" if out_file.endswith("gz") else ""
            sample_name_str = items[0]["name"][1]
            cl2 = ("cat {scalpel_tmp_file} | sed 's/sample_name/{sample_name_str}/g' | {vcfallelicprimitives} | {vcfstreamsort} {compress_cmd} > {tx_out_file}")
            do.run(cl2.format(**locals()), "Finalising Scalpel variants", {})
    ann_file = annotation.annotate_nongatk_vcf(out_file, align_bams,
                                               assoc_files["dbsnp"],
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
        with file_transaction(out_file) as tx_out_file:
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
            min_cov = "5" # minimum coverage (default 5)
            opts += " --mincov %s" % min_cov
            cl = ("{scalpel} --somatic {opts} --tumor {paired.tumor_bam} --normal {paired.normal_bam}")
            bam.index(paired.tumor_bam, config)
            bam.index(paired.normal_bam, config)
            do.run(cl.format(**locals()), "Genotyping paired variants with Scalpel", {})
            scalpel_tmp_file = os.path.join(tmp_path, "main/somatic."+min_cov+"x.indel.vcf")
            scalpel_tmp_file_common = os.path.join(tmp_path, "main/common."+min_cov+"x.indel.vcf")
            compress_cmd = "| bgzip -c" if out_file.endswith("gz") else ""
            cl2 = ("cat {scalpel_tmp_file} <(grep -vE '^#' {scalpel_tmp_file_common} | sed 's/PASS/REJECT/g') | sed 's/sample_name/{paired.tumor_name}/g' | {vcfstreamsort} {compress_cmd} > {tx_out_file}")
            do.run(cl2.format(**locals()), "Finalising Scalpel variants", {})
    ann_file = annotation.annotate_nongatk_vcf(out_file, align_bams,
                                               assoc_files["dbsnp"], ref_file,
                                               config)
    return ann_file