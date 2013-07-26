"""Provide variant calling with VarScan from TGI at Wash U.

http://varscan.sourceforge.net/
"""
import contextlib
import itertools
import os
import tempfile

from bcbio.pipeline import config_utils
from bcbio.provenance import do, programs
from bcbio.variation import samtools
from bcbio.variation.genotype import write_empty_vcf
from bcbio.variation.vcfutils import combine_variant_files

import pysam


def run_varscan(align_bams, items, ref_file, assoc_files,
                region=None, out_file=None):
    call_file = samtools.shared_variantcall(_varscan_work, "varscan", align_bams,
                                            ref_file, items[0]["config"], assoc_files, region, out_file)
    return call_file


def run_varscan_paired(align_bams, items, ref_file, assoc_files,
                       region=None, out_file=None):

    call_file = samtools.shared_variantcall(_varscan_paired, "varscan",
                                            align_bams, items[0]["config"],
                                            assoc_files, region, out_file)

    return call_file


def _varscan_paired(align_bams, items, ref_file, target_regions, out_file):

    max_read_depth = "1000"
    config = items[0]["config"]

    version = programs.jar_versioner("varscan", "VarScan")(config)
    if version < "v2.3.5":
        raise IOError(
            "Please install version 2.3.5 or better of VarScan with support "
            "for multisample calling and indels in VCF format.")
    varscan_jar = config_utils.get_jar(
        "VarScan",
        config_utils.get_program("varscan", config, "dir"))

    resources = config_utils.get_resources("varscan", config)
    jvm_opts = " ".join(resources.get("jvm_opts", ["-Xmx750m", "-Xmx2g"]))
    remove_zerocoverage = "grep -v -P '\t0\t\t$'"

    with (tempfile.NamedTemporaryFile(), tempfile.NamedTemporaryFile()) as (
        normal_tmp_mpileup, tumor_tmp_mpilpeup):
        for bamfile, item in itertools.izip(align_bams, items):

            metadata = item["metadata"]

            mpileup = samtools.prep_mpileup([bamfile], ref_file,
                                            max_read_depth, config,
                                            target_regions=target_regions,
                                            want_bcf=False)

            if metadata["phenotype"] == "normal":
                cmd = "{mpileup} | {remove_zerocoverage} > {normal_tmp_pileup}"
            elif metadata["phenotype"] == "tumor":
                cmd = "{mpileup} | {remove_zerocoverage} > {tumor_tmp_pileup}"

            cmd = cmd.format(**locals())
            do.run(cmd)
            # FIXME How to check for success?

        varscan_cmd = ("java {jvm_opts} -jar {varscan_jar} somatic"
                       " {normal_tmp_mpileup} {tumor_tmp_mpileup} {out_file}"
                       "--output-vcf --min-coverage 5 --p-value 0.98")
        #FIXME This currently generates *two* files that need to be merged!

        out_file_snp = out_file + ".snp"
        out_file_indel = out_file + ".indel"

        out_file = combine_variant_files([out_file_snp, out_file_indel],
                                         out_file, ref_file, config,
                                         region=target_regions)

        os.remove(out_file_snp)
        os.remove(out_file_indel)

        if os.path.getsize(out_file) == 0:
            write_empty_vcf(out_file)


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
        write_empty_vcf(out_file)
