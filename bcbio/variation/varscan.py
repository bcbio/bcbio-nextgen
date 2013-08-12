"""Provide variant calling with VarScan from TGI at Wash U.

http://varscan.sourceforge.net/
"""
import contextlib
import os

from bcbio.pipeline import config_utils
from bcbio.provenance import do, programs
from bcbio.variation import samtools, vcfutils

import pysam

def run_varscan(align_bams, items, ref_file, assoc_files,
                region=None, out_file=None):
    call_file = samtools.shared_variantcall(_varscan_work, "varscan", align_bams,
                                            ref_file, items[0]["config"], assoc_files, region, out_file)
    return call_file

def _create_sample_list(in_bams, vcf_file):
    """Pull sample names from input BAMs and create input sample list.
    """
    out_file = "%s-sample_list.txt" % os.path.splitext(vcf_file)[0]
    with open(out_file, "w") as out_handle:
        for in_bam in in_bams:
            with contextlib.closing(pysam.Samfile(in_bam, "rb")) as work_bam:
                for rg in work_bam.header.get("RG", []):
                    out_handle.write("%s\n" % rg["SM"])
    return out_file

def _varscan_work(align_bams, ref_file, config, target_regions, out_file):
    """Perform SNP and indel genotyping with VarScan.
    """
    max_read_depth = "1000"
    version = programs.jar_versioner("varscan", "VarScan")(config)
    if version < "v2.3.5":
        raise IOError("Please install version 2.3.5 or better of VarScan with support "
                      "for multisample calling and indels in VCF format.")
    varscan_jar = config_utils.get_jar("VarScan",
                                       config_utils.get_program("varscan", config, "dir"))
    resources = config_utils.get_resources("varscan", config)
    jvm_opts = " ".join(resources.get("jvm_opts", ["-Xmx750m", "-Xmx2g"]))
    sample_list = _create_sample_list(align_bams, out_file)
    mpileup = samtools.prep_mpileup(align_bams, ref_file, max_read_depth, config,
                                    target_regions=target_regions, want_bcf=False)
    # VarScan fails to generate a header on files that start with
    # zerocoverage calls; strip these with grep, we're not going to
    # call on them
    remove_zerocoverage = "grep -v -P '\t0\t\t$'"
    cmd = ("{mpileup} | {remove_zerocoverage} "
           "| java {jvm_opts} -jar {varscan_jar} mpileup2cns --min-coverage 5 --p-value 0.98 "
           "  --vcf-sample-list {sample_list} --output-vcf --variants "
           "> {out_file}")
    cmd = cmd.format(**locals())
    do.run(cmd, "Varscan".format(**locals()), None,
           [do.file_exists(out_file)])
    os.remove(sample_list)
    # VarScan can create completely empty files in regions without
    # variants, so we create a correctly formatted empty file
    if os.path.getsize(out_file) == 0:
        vcfutils.write_empty_vcf(out_file)
