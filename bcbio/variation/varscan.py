"""Provide variant calling with VarScan from TGI at Wash U.

http://varscan.sourceforge.net/
"""
import contextlib
import itertools
import os

from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.provenance import do, programs
from bcbio.utils import file_exists
from bcbio.variation import samtools
from bcbio.variation.vcfutils import (combine_variant_files, write_empty_vcf,
                                      is_sample_pair, get_paired_bams)

import pysam


def run_varscan(align_bams, items, ref_file, assoc_files,
                region=None, out_file=None):

    if len(align_bams) == 2 and all(item["metadata"].get("phenotype")
                                    is not None for item in items):
        call_file = samtools.shared_variantcall(_varscan_paired, "varscan",
                                            align_bams, ref_file, items,
                                            assoc_files, region, out_file)
    else:
        call_file = samtools.shared_variantcall(_varscan_work, "varscan",
                                                align_bams, ref_file,
                                                items, assoc_files,
                                                region, out_file)
    return call_file


def _varscan_paired(align_bams, ref_file, items, target_regions, out_file):

    """Run a paired VarScan analysis, also known as "somatic". """

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

    # No need for names in VarScan, hence the "_"

    tumor_bam, _, normal_bam, _ = get_paired_bams(align_bams, items)

    if not file_exists(out_file):
        base, ext = os.path.splitext(out_file)
        cleanup_files = []
        for fname, mpext in [(normal_bam, "normal"), (tumor_bam, "tumor")]:
            mpfile = "%s-%s.mpileup" % (base, mpext)
            cleanup_files.append(mpfile)
            with file_transaction(mpfile) as mpfile_tx:
                mpileup = samtools.prep_mpileup([mpfile_tx], ref_file,
                                                max_read_depth, config,
                                                target_regions=target_regions,
                                                want_bcf=False)
                cmd = "{mpileup} | {remove_zerocoverage} > {mpfile_tx}"
                cmd = cmd.format(**locals())
                do.run(cmd)

        # First index is normal, second is tumor
        normal_tmp_mpileup = cleanup_files[0]
        tumor_tmp_mpileup = cleanup_files[1]

        varscan_cmd = ("java {jvm_opts} -jar {varscan_jar} somatic"
                       " {normal_tmp_mpileup} {tumor_tmp_mpileup} {base}"
                       "--output-vcf --min-coverage 5 --p-value 0.98")

        indel_file = base + ".indel"
        snp_file = base + ".snp"

        cleanup_files.append(indel_file)
        cleanup_files.append(snp_file)

        with (file_transaction(indel_file), file_transaction(snp_file)) as (
            tx_indel, tx_snp):
            varscan_cmd = varscan_cmd.format(**locals())
            do.run(varscan_cmd)

        out_file = combine_variant_files([snp_file, indel_file],
                                         out_file, ref_file, config,
                                         region=target_regions)

        # Remove cleanup files

        for extra_file in cleanup_files:
            os.remove(extra_file)

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


def _varscan_work(align_bams, ref_file, items, target_regions, out_file):
    """Perform SNP and indel genotyping with VarScan.
    """

    config = items[0]["config"]

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
