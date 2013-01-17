"""Provide variant calling with VarScan from TGI at Wash U.

http://varscan.sourceforge.net/
"""
import os
import shutil
import contextlib

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.variation import samtools

import sh
import pysam

def run_varscan(align_bams, ref_file, config,
                dbsnp=None, region=None, out_file=None):
    call_file = samtools.shared_variantcall(_varscan_work, "varscan", align_bams,
                                            ref_file, config, dbsnp, region, out_file)
    _fix_varscan_vcf(call_file, align_bams)
    return call_file

def _fix_varscan_vcf(orig_file, in_bams):
    """Fixes issues with the standard VarScan VCF output.

    - Remap sample names back to those defined in the input BAM file.
    - Convert indels into correct VCF representation.
    """
    tmp_file = utils.append_stem(orig_file, "origsample", "-")
    if not utils.file_exists(tmp_file):
        shutil.move(orig_file, tmp_file)
        with file_transaction(orig_file) as tx_out_file:
            with open(tmp_file) as in_handle:
                with open(orig_file, "w") as out_handle:
                    for line in in_handle:
                        parts = line.split("\t")
                        if line.startswith("#CHROM"):
                            line = _fix_sample_line(line, in_bams)
                        elif not line.startswith("#") and parts[4].startswith(("+", "-")):
                            line = _fix_indel_line(parts)
                        out_handle.write(line)

def _fix_indel_line(parts):
    """Convert VarScan indel representations into standard VCF.
    """
    ref = parts[3]
    alt = parts[4]
    mod_alt = alt[0]
    seq_alt = alt[1:]
    if mod_alt == "+":
        new_ref = ref
        new_alt = ref + seq_alt
    elif mod_alt == "-":
        new_ref = ref + seq_alt
        new_alt = ref
    parts[3] = new_ref
    parts[4] = new_alt
    return "\t".join(parts)

def _fix_sample_line(line, in_bams):
    """Pull sample names from input BAMs and replace VCF file header.
    """
    samples = []
    for in_bam in in_bams:
        with contextlib.closing(pysam.Samfile(in_bam, "rb")) as work_bam:
            for rg in work_bam.header.get("RG", []):
                samples.append(rg["SM"])
    parts = line.split("\t")
    standard = parts[:9]
    old_samples = parts[9:]
    if len(old_samples) == 0:
        return line
    else:
        assert len(old_samples) == len(samples), (old_samples, samples)
        return "\t".join(standard + samples) + "\n"

def _varscan_work(align_bams, ref_file, config, target_regions, out_file):
    """Perform SNP and indel genotyping with VarScan.
    """
    max_read_depth = 1000
    varscan_jar = config_utils.get_jar("VarScan",
                                       config_utils.get_program("varscan", config, "dir"))
    with open(out_file, "w") as out_handle:
        mpileup = samtools.prep_mpileup(align_bams, ref_file, max_read_depth, target_regions,
                                        want_bcf=False)
        varscan = sh.Command("java").bake("-jar", varscan_jar, "mpileup2cns",
                                          "--min-coverage", "5",
                                          "--p-value", "0.98",
                                          "--output-vcf", "--variants", _out=out_handle)
        varscan(mpileup())
