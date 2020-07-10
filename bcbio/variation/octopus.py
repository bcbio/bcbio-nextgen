"""Octopus germline and somatic calling.

https://github.com/luntergroup/octopus
"""
import os
import subprocess

import pysam
import six

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import shared
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import bamprep, bedutils, vcfutils

def run(align_bams, items, ref_file, assoc_files, region, out_file):
    """Run octopus variant calling, handling both somatic and germline calling.
    """
    if not utils.file_exists(out_file):
        paired = vcfutils.get_paired_bams(align_bams, items)
        regions = _get_regions(region, out_file, items)
        if paired:
            return _run_somatic(paired, ref_file, regions, out_file)
        else:
            return _run_germline(align_bams, items, ref_file, regions, out_file)
    return out_file

def _get_regions(region, out_file, items):
    """Retrieve region to run analysis in. Handles no targets, BED and regions.
    """
    vrs = bedutils.population_variant_regions(items)
    target = shared.subset_variant_regions(vrs, region, out_file, items=items, do_merge=True)
    if target:
        if isinstance(target, six.string_types) and os.path.isfile(target):
            return "--regions-file %s" % target
        else:
            return "--regions %s" % bamprep.region_to_gatk(target)
    else:
        return ""

def _produce_compatible_vcf(out_file, data, is_somatic=False):
    """Create a compatible VCF that downstream tools can deal with.

    - htsjdk and thus GATK and Picard do not support VCF4.3:
      https://github.com/broadinstitute/gatk/issues/2092
    - Use octopus legacy format to avoid incompatibilities.
      https://github.com/luntergroup/octopus#output-format
    - Fixes `##contig` lines since octopus only writes contigs
      used in the BED file region, causing incompatibilies with
      GatherVcfs when merging
    - Fixes alleles prefixed with '*' like 'C,*T' which cause
      downstream failures when reading with GATK.
    """
    base, ext = utils.splitext_plus(out_file)
    legacy_file = "%s.legacy%s" % (base, ext)
    if is_somatic:
        legacy_file = _covert_to_diploid(legacy_file, data)
    final_file = "%s.vcf.gz" % base
    cat_cmd = "zcat" if legacy_file.endswith(".gz") else "cat"
    contig_cl = vcfutils.add_contig_to_header_cl(dd.get_ref_file(data), out_file)
    remove_problem_alleles = r"sed 's/,\*\([A-Z]\)/,\1/'"
    cmd = ("{cat_cmd} {legacy_file} | sed 's/fileformat=VCFv4.3/fileformat=VCFv4.2/' | "
           "{remove_problem_alleles} | {contig_cl} | bgzip -c > {final_file}")
    do.run(cmd.format(**locals()), "Produce compatible VCF output file from octopus")
    return vcfutils.bgzip_and_index(out_file, data["config"])

def _covert_to_diploid(in_file, data):
    """Converts non-diploid somatic outputs into diploid.

    https://github.com/luntergroup/octopus/wiki/Case-study:-Tumour-only-UMI#evaluate-variant-calls
    """
    sample = dd.get_sample_name(data)
    out_file = "%s-diploid.vcf" % utils.splitext_plus(in_file)[0]
    in_vcf = pysam.VariantFile(in_file)
    out_vcf = pysam.VariantFile(out_file, 'w', header=in_vcf.header)
    for record in in_vcf:
        gt = list(record.samples[sample]['GT'])
        if 'SOMATIC' in record.info:
            for allele in set(gt):
                if allele != gt[0]:
                    record.samples[sample]['GT'] = gt[0], allele
                    out_vcf.write(record)
        else:
            if len(gt) == 1:
                record.samples[sample]['GT'] = gt
            else:
                record.samples[sample]['GT'] = gt[0], gt[1]
            out_vcf.write(record)
    in_vcf.close()
    out_vcf.close()
    return vcfutils.bgzip_and_index(out_file, data["config"])

def _run_germline(align_bams, items, ref_file, regions, out_file):
    """Run germline calling, handling populations.

    TODO:
      - We could better handle trio calling with ped inputs as octopus
        has special support.
    """
    align_bams = " ".join(align_bams)
    cores = dd.get_num_cores(items[0])
    cmd = ("octopus --threads {cores} --reference {ref_file} --reads {align_bams} "
           "{regions} "
           "--working-directory {tmp_dir} "
           "-o {tx_out_file} --legacy")
    with file_transaction(items[0], out_file) as tx_out_file:
        tmp_dir = os.path.dirname(tx_out_file)
        do.run(cmd.format(**locals()), "Octopus germline calling")
        _produce_compatible_vcf(tx_out_file, items[0])
    return out_file

def _run_somatic(paired, ref_file, regions, out_file):
    """Run somatic calling with octopus, handling both paired and tumor-only cases.

    Tweaks for low frequency, tumor only and UMI calling documented in:
    https://github.com/luntergroup/octopus/blob/develop/configs/UMI.config
    """
    align_bams = paired.tumor_bam
    if paired.normal_bam:
        align_bams += " %s --normal-sample %s" % (paired.normal_bam, paired.normal_name)
    cores = dd.get_num_cores(paired.tumor_data)
    # Do not try to search below 0.4% currently as leads to long runtimes
    # https://github.com/luntergroup/octopus/issues/29#issuecomment-428167979
    min_af = max([float(dd.get_min_allele_fraction(paired.tumor_data)) / 100.0, 0.004])
    min_af_floor = min_af / 4.0
    cmd = ("octopus --threads {cores} --reference {ref_file} --reads {align_bams} "
           "{regions} "
           "--min-credible-somatic-frequency {min_af_floor} --min-expected-somatic-frequency {min_af} "
           "--downsample-above 4000 --downsample-target 4000 --min-kmer-prune 5 --min-bubble-score 20 "
           "--max-haplotypes 200 --somatic-snv-mutation-rate '5e-4' --somatic-indel-mutation-rate '1e-05' "
           "--target-working-memory 5G --target-read-buffer-footprint 5G --max-somatic-haplotypes 3 "
           "--caller cancer "
           "--working-directory {tmp_dir} "
           "-o {tx_out_file} --legacy")
    if not paired.normal_bam:
        cmd += (" --tumour-germline-concentration 5")
    if dd.get_umi_type(paired.tumor_data) or _is_umi_consensus_bam(paired.tumor_bam):
        cmd += (" --allow-octopus-duplicates --overlap-masking 0 "
                "--somatic-filter-expression 'GQ < 200 | MQ < 30 | SB > 0.2 | SD[.25] > 0.1 "
                "| BQ < 40 | DP < 100 | MF > 0.1 | AD < 5 | CC > 1.1 | GQD > 2'")
    with file_transaction(paired.tumor_data, out_file) as tx_out_file:
        tmp_dir = os.path.dirname(tx_out_file)
        do.run(cmd.format(**locals()), "Octopus somatic calling")
        _produce_compatible_vcf(tx_out_file, paired.tumor_data, is_somatic=True)
    return out_file

def _is_umi_consensus_bam(in_file):
    """Check if input BAM file generated by fgbio consensus calls on UMIs.

    Identify these by lack of duplicated reads.
    This is useful for pre-aligned consensus BAMs feeding into octopus.
    """
    cmd = "samtools view -h %s | head -500000 | samtools view -c -f 1024"
    count = subprocess.check_output(cmd % in_file, shell=True)
    return int(count) == 0
