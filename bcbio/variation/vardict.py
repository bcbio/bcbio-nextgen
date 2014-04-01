"""Sensitive variant calling using VarDict.

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


def region_to_vardict(region):
    if isinstance(region, (list, tuple)):
        chrom, start, end = region
        return "%s:%s-%s" % (chrom, start, end)
    else:
        return region

def _vardict_options_from_config(items, config, out_file, region=None):
    opts = []

    resources = config_utils.get_resources("vardict", config)
    if resources.get("options"):
        opts += resources["options"]

    variant_regions = utils.get_in(config, ("algorithm", "variant_regions"))
    target = subset_variant_regions(variant_regions, region, out_file)
    if target:
        if isinstance(target, basestring) and os.path.isfile(target):
            opts += [target] # this must be the last option
        else:
            opts += ["-R", region_to_vardict(target)]
    return opts

def run_vardict(align_bams, items, ref_file, assoc_files, region=None,
                  out_file=None):
    """Run VarDict variant calling.
    """
    call_file = _run_vardict_caller(align_bams, items, ref_file,
                                          assoc_files, region, out_file)

    return call_file

def _run_vardict_caller(align_bams, items, ref_file, assoc_files,
                          region=None, out_file=None):
    """Detect SNPs and indels with VarDict.
    """
    config = items[0]["config"]
    if out_file is None:
        out_file = "%s-variants.vcf.gz" % os.path.splitext(align_bams[0])[0]
    if not utils.file_exists(out_file):
        tmp_out1 = out_file + ".temp1.txt"
        tmp_out2 = out_file + ".temp2.txt"
        with file_transaction(out_file) as tx_out_file:
            for align_bam in align_bams:
                bam.index(align_bam, config)
            
            #vardict = config_utils.get_program("vardict", config)
            input_bams = " ".join("-b %s" % x for x in align_bams)
            
            opts = " ".join(_vardict_options_from_config(items, config, out_file, region))            
            vcfstreamsort = config_utils.get_program("vcfstreamsort", config)
            compress_cmd = "| bgzip -c" if out_file.endswith("gz") else ""            
            freq = 0.01
            item = 
            sample = item["name"][1]
            input_bam = 
            if "GRC" in ...:
            with file_transaction(tmp_out1) as tx_tmp1_file:
                cmd = ("vardict.pl -c 1 -S 2 -E 3 -g 4 -x 0 -G {ref_file} -f {freq} -k 3 " 
                       "-X 5 -N {sample} -b {input_bam} {opts} > {tx_tmp1_file}")
                do.run(cmd.format(**locals()), "Genotyping with VarDict: Inference", {})
            with file_transaction(tmp_out2) as tx_tmp2_file:
                cmd = ("teststrandbias.R {tmp_out} > {tx_tmp2_file}")
                do.run(cmd.format(**locals()), "Genotyping with VarDict: Strand bias", {})
            cmd = ("var2vcf_valid.pl -S -f {freq} {tmp_out2} | {vcfstreamsort} | {compress_cmd} > {tx_out_file}")
            do.run(cmd.format(**locals()), "Genotyping with VarDict: VCF output", {})
    ann_file = annotation.annotate_nongatk_vcf(out_file, align_bams,
                                               assoc_files["dbsnp"],
                                               ref_file, config)
    return ann_file

