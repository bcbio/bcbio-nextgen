"""Compare two input files, with optional BED exclusion regions using SMaSH calldiff:

https://github.com/kwestbrooks/smash/tree/master/calldiff

Currently a work in progress -- not yet integrated.

Testing Usage:
  compare_calldiff.py <ref_file> <truth_vcf> <eval_vcf> [<one or more BED files]
"""
import os
import sys

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.variation import vcfutils

def compare(ref_file, truth_vcf, eval_vcf, *bed_files):
    out_dir = utils.safe_makedir(os.path.join(os.getcwd(), "validate"))
    region_bed = _final_bed_region(bed_files, eval_vcf, out_dir)
    if region_bed:
        truth_vcf = _subset_vcf(truth_vcf, region_bed, out_dir)
        eval_vcf = _subset_vcf(eval_vcf, region_bed, out_dir)
    _do_smash_calldiff(truth_vcf, eval_vcf, ref_file, out_dir)

def _do_smash_calldiff(truth_vcf, eval_vcf, ref_file, out_dir):
    cmd = ["java", "-jar", "calldiff-jar-with-dependencies.jar",
           "--lhs_vcf", truth_vcf, "--rhs_vcf", eval_vcf,
           "--reference_fasta", ref_file, "--presorted"]
    do.run(cmd, "Compare files with SMaSH calldiff")

def _subset_vcf(vcf_file, bed_file, out_dir):
    """Restrict VCF to only regions defined in the initial input BED file.
    """
    #out_file = os.path.join(out_dir, "%s-cmp.vcf.gz" % utils.splitext_plus(os.path.basename(vcf_file))[0])
    out_file = os.path.join(out_dir, "%s-cmp.vcf" % utils.splitext_plus(os.path.basename(vcf_file))[0])
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            #cmd = ("bedtools intersect -a {vcf_file} -b {bed_file} -wa -header | "
            #       "bgzip -c > {tx_out_file}")
            cmd = ("bedtools intersect -a {vcf_file} -b {bed_file} -wa -header > {tx_out_file}")
            do.run(cmd.format(**locals()), "Subset inputs to BED file of interest")
    #return vcfutils.bgzip_and_index(out_file, {})
    return out_file

def _final_bed_region(bed_files, eval_vcf, out_dir):
    """Prepare final BED region: combined intersection of all input BEDs.
    """
    if len(bed_files) == 0:
        return None
    elif len(bed_files) == 1:
        return bed_files[0]
    else:
        out_file = os.path.join(out_dir,
                                "%s-regions.bed.gz" % utils.splitext_plus(os.path.basename(eval_vcf))[0])
        if not utils.file_exists(out_file):
            with file_transaction(out_file) as tx_out_file:
                if len(bed_files) == 2:
                    cmd = ("bedtools intersect -a {bed_files[0]} -b {bed_files[1]} | "
                           "sort -k1,1 -k2,2n | bgzip -c > {tx_out_file}")
                else:
                    cmd = ("bedtools multiinter -i {" ".join(bed_files)} | "
                           "sort -k1,1 -k2,2n | bgzip -c > {tx_out_file}")
                do.run(cmd.format(**locals()), "Prepare combined BED region")
        return vcfutils.bgzip_and_index(out_file, {})

if __name__ == "__main__":
    compare(*sys.argv[1:])
