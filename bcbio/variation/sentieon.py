"""Variant calling support for Sentieon tools.

Sentieon provides optimized versions of standard tools like GATK HaplotypeCaller
and MuTect2 as well as their own developed versions. These require a license
from Sentieon for use:

http://sentieon.com/about/
https://peerj.com/preprints/1672/
"""
import os
import pprint

import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils, shared
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import bamprep, bedutils, joint, vcfutils

import six


def license_export(data):
    """Retrieve export statement for sentieon license server.
    """
    resources = config_utils.get_resources("sentieon", data["config"])
    server = resources.get("keyfile")
    if not server:
        server = tz.get_in(["resources", "sentieon", "keyfile"], data)
    if not server:
        raise ValueError("Need to set resources keyfile with URL:port of license server, local license file or "
                         "environmental variables to export \n"
                         "http://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html#resources\n"
                         "Configuration: %s" % pprint.pformat(data))
    if isinstance(server, six.string_types):
        return "export SENTIEON_LICENSE=%s && " % server
    else:
        assert isinstance(server, dict), server
        exports = ""
        for key, val in server.items():
            exports += "export %s=%s && " % (key.upper(), val)
        return exports

def _get_interval(variant_regions, region, out_file, items):
    """Retrieve interval to run analysis in. Handles no targets, BED and regions

    region can be a single region or list of multiple regions for multicore calling.
    """
    target = shared.subset_variant_regions(variant_regions, region, out_file, items)
    if target:
        if isinstance(target, six.string_types) and os.path.isfile(target):
            return "--interval %s" % target
        else:
            return "--interval %s" % bamprep.region_to_gatk(target)
    else:
        return ""

def run_tnscope(align_bams, items, ref_file, assoc_files,
                     region=None, out_file=None):
    """Call variants with Sentieon's TNscope somatic caller.
    """
    if out_file is None:
        out_file = "%s-variants.vcf.gz" % utils.splitext_plus(align_bams[0])[0]
    if not utils.file_exists(out_file):
        variant_regions = bedutils.population_variant_regions(items, merged=True)
        interval = _get_interval(variant_regions, region, out_file, items)
        with file_transaction(items[0], out_file) as tx_out_file:
            paired = vcfutils.get_paired_bams(align_bams, items)
            assert paired and paired.normal_bam, "Require normal BAM for Sentieon TNscope"
            dbsnp = "--dbsnp %s" % (assoc_files.get("dbsnp")) if "dbsnp" in assoc_files else ""
            license = license_export(items[0])
            cores = dd.get_num_cores(items[0])
            cmd = ("{license}sentieon driver -t {cores} -r {ref_file} "
                   "-i {paired.tumor_bam} -i {paired.normal_bam} {interval} "
                   "--algo TNscope "
                   "--tumor_sample {paired.tumor_name} --normal_sample {paired.normal_name} "
                   "{dbsnp} {tx_out_file}")
            do.run(cmd.format(**locals()), "Sentieon TNscope")
    return out_file

def run_tnhaplotyper(align_bams, items, ref_file, assoc_files,
                     region=None, out_file=None):
    """Call variants with Sentieon's TNhaplotyper (MuTect2 like).
    """
    if out_file is None:
        out_file = "%s-variants.vcf.gz" % utils.splitext_plus(align_bams[0])[0]
    if not utils.file_exists(out_file):
        variant_regions = bedutils.population_variant_regions(items, merged=True)
        interval = _get_interval(variant_regions, region, out_file, items)
        with file_transaction(items[0], out_file) as tx_out_file:
            paired = vcfutils.get_paired_bams(align_bams, items)
            assert paired.normal_bam, "Require normal BAM for Sentieon TNhaplotyper"
            dbsnp = "--dbsnp %s" % (assoc_files.get("dbsnp")) if "dbsnp" in assoc_files else ""
            cosmic = "--cosmic %s" % (assoc_files.get("cosmic")) if "cosmic" in assoc_files else ""
            license = license_export(items[0])
            tx_orig_file = "%s-orig%s" % utils.splitext_plus(tx_out_file)
            cores = dd.get_num_cores(items[0])
            cmd = ("{license}sentieon driver -t {cores} -r {ref_file} "
                   "-i {paired.tumor_bam} -i {paired.normal_bam} {interval} "
                   "--algo TNhaplotyper "
                   "--tumor_sample {paired.tumor_name} --normal_sample {paired.normal_name} "
                   "{dbsnp} {cosmic} {tx_orig_file}")
            do.run(cmd.format(**locals()), "Sentieon TNhaplotyper")
            cmd = ("gunzip -c {tx_orig_file} | "
                   "sed 's/ID=ECNT,Number=1,Type=Integer/ID=ECNT,Number=1,Type=String/' | "
                   "sed 's/ID=HCNT,Number=1,Type=Integer/ID=HCNT,Number=1,Type=String/' | "
                   "sed 's/ID=NLOD,Number=1,Type=Float/ID=NLOD,Number=1,Type=String/' | "
                   "sed 's/ID=TLOD,Number=1,Type=Float/ID=TLOD,Number=1,Type=String/' | "
                   "sed 's/ID=PON,Number=1,Type=Integer/ID=PON,Number=1,Type=String/' | "
                   "bgzip -c > {tx_out_file}")
            do.run(cmd.format(**locals()), "Sentieon TNhaplotyper: make headers GATK compatible")
            vcfutils.bgzip_and_index(tx_out_file, items[0]["config"])
    return out_file

def run_haplotyper(align_bams, items, ref_file, assoc_files,
                     region=None, out_file=None):
    """Call variants with Sentieon's haplotyper (GATK HaplotypeCaller like).
    """
    if out_file is None:
        out_file = "%s-variants.vcf.gz" % utils.splitext_plus(align_bams[0])[0]
    if not utils.file_exists(out_file):
        variant_regions = bedutils.population_variant_regions(items, merged=True)
        interval = _get_interval(variant_regions, region, out_file, items)
        with file_transaction(items[0], out_file) as tx_out_file:
            dbsnp = "--dbsnp %s" % (assoc_files.get("dbsnp")) if "dbsnp" in assoc_files else ""
            bams = " ".join(["-i %s" % x for x in align_bams])
            license = license_export(items[0])
            cores = dd.get_num_cores(items[0])
            out_mode = "--emit_mode gvcf" if joint.want_gvcf(items) else ""
            cmd = ("{license}sentieon driver -t {cores} -r {ref_file} "
                   "{bams} {interval} --algo Haplotyper {out_mode} {dbsnp} {tx_out_file}")
            do.run(cmd.format(**locals()), "Sentieon Haplotyper")
    return out_file

def run_gvcftyper(vrn_files, out_file, region, data):
    """Produce joint called variants from input gVCF files.
    """
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            license = license_export(data)
            ref_file = dd.get_ref_file(data)
            input_files = " ".join(vrn_files)
            region = bamprep.region_to_gatk(region)
            cmd = ("{license}sentieon driver -r {ref_file} --interval {region} "
                   "--algo GVCFtyper {tx_out_file} {input_files}")
            do.run(cmd.format(**locals()), "Sentieon GVCFtyper")
    return out_file

def bqsr_table(data):
    """Generate recalibration tables as inputs to BQSR.
    """
    in_file = dd.get_align_bam(data)
    out_file = "%s-recal-table.txt" % utils.splitext_plus(in_file)[0]
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(data, out_file) as tx_out_file:
            assoc_files = dd.get_variation_resources(data)
            known = "-k %s" % (assoc_files.get("dbsnp")) if "dbsnp" in assoc_files else ""
            license = license_export(data)
            cores = dd.get_num_cores(data)
            ref_file = dd.get_ref_file(data)
            cmd = ("{license}sentieon driver -t {cores} -r {ref_file} "
                   "-i {in_file} --algo QualCal {known} {tx_out_file}")
            do.run(cmd.format(**locals()), "Sentieon QualCal generate table")
    return out_file

def apply_bqsr(data):
    """Apply recalibration, producing a updated BAM file.
    """
    in_file = dd.get_align_bam(data)
    out_table_file = "%s-recal-table-post.txt" % utils.splitext_plus(in_file)[0]
    out_file = "%s-recal.bam" % utils.splitext_plus(in_file)[0]
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(data, out_file, out_table_file) as (tx_out_file, tx_table_file):
            assoc_files = dd.get_variation_resources(data)
            known = "-k %s" % (assoc_files.get("dbsnp")) if "dbsnp" in assoc_files else ""
            license = license_export(data)
            cores = dd.get_num_cores(data)
            ref_file = dd.get_ref_file(data)
            cmd = ("{license}sentieon driver -t {cores} -r {ref_file} "
                   "-i {in_file} --algo QualCal {known} {tx_table_file} "
                   "--algo ReadWriter {tx_out_file}")
            do.run(cmd.format(**locals()), "Sentieon QualCal apply recalibration")
    return out_file
