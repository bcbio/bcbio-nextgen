"""Octopus germline and somatic calling.

https://github.com/luntergroup/octopus
"""
import os

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import shared
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import bedutils, vcfutils

def run(align_bams, items, ref_file, assoc_files, region, out_file):
    """Run octopus variant calling, handling both somatic and germline calling.
    """
    if not utils.file_exists(out_file):
        paired = vcfutils.get_paired_bams(align_bams, items)
        vrs = bedutils.population_variant_regions(items)
        target = shared.subset_variant_regions(vrs, region,
                                               out_file, items=items, do_merge=True)
        if paired:
            return _run_somatic(paired, ref_file, target, out_file)
        else:
            return _run_germline(align_bams, items, ref_file, target, out_file)
    return out_file

def _produce_compatible_vcf(out_file, data):
    """Create a compatible VCF that downstream tools can deal with.

    - htsjdk and thus GATK and Picard do not support VCF4.3:
      https://github.com/broadinstitute/gatk/issues/2092
    - Use octopus legacy format to avoid incompatibilities.
      https://github.com/luntergroup/octopus#output-format
    - Fixes `##contig` lines since octopus only writes contigs
      used in the BED file region, causing incompatibilies with
      GatherVcfs when merging
    """
    base, ext = utils.splitext_plus(out_file)
    legacy_file = "%s.legacy%s" % (base, ext)
    final_file = "%s.vcf.gz" % base
    cat_cmd = "zcat" if legacy_file.endswith(".gz") else "cat"
    contig_cl = vcfutils.add_contig_to_header_cl(dd.get_ref_file(data), out_file)
    cmd = ("{cat_cmd} {legacy_file} | sed 's/fileformat=VCFv4.3/fileformat=VCFv4.2/' | "
           "{contig_cl} | bgzip -c > {final_file}")
    do.run(cmd.format(**locals()), "Produce compatible VCF output file from octopus")
    return vcfutils.bgzip_and_index(out_file, data["config"])

def _run_germline(align_bams, items, ref_file, target, out_file):
    """Run germline calling, handling populations.

    TODO:
      - We could better handle trio calling with ped inputs as octopus
        has special support.
    """
    align_bams = " ".join(align_bams)
    cores = dd.get_num_cores(items[0])
    cmd = ("octopus --threads {cores} --reference {ref_file} --reads {align_bams} "
           "--regions-file {target} "
           "--working-directory {tmp_dir} "
           "-o {tx_out_file} --legacy")
    with file_transaction(items[0], out_file) as tx_out_file:
        tmp_dir = os.path.dirname(tx_out_file)
        do.run(cmd.format(**locals()), "Octopus germline calling")
        _produce_compatible_vcf(tx_out_file, items[0])
    return out_file

def _run_somatic(paired, ref_file, target, out_file):
    """Run somatic calling with octopus, handling both paired and tumor-only cases.

    TODO:
      - Should we set downsample-above and downsample-target which default to 1000
        and 500? How will these effect high depth panels and runtimes?
    """
    align_bams = paired.tumor_bam
    if paired.normal_bam:
        align_bams += " %s --normal-sample %s" % (paired.normal_bam, paired.normal_name)
    cores = dd.get_num_cores(paired.tumor_data)
    min_af = float(dd.get_min_allele_fraction(paired.tumor_data)) / 100.0
    cmd = ("octopus --threads {cores} --reference {ref_file} --reads {align_bams} "
           "--regions-file {target} "
           "--min-credible-somatic-frequency {min_af} "
           "-C cancer "
           "--working-directory {tmp_dir} "
           "-o {tx_out_file} --legacy")
    with file_transaction(paired.tumor_data, out_file) as tx_out_file:
        tmp_dir = os.path.dirname(tx_out_file)
        do.run(cmd.format(**locals()), "Octopus somatic calling")
        _produce_compatible_vcf(tx_out_file, paired.tumor_data)
    return out_file
