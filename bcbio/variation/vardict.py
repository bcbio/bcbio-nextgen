"""Sensitive variant calling using VarDict.

Defaults to using the faster, equally sensitive Java port:

https://github.com/AstraZeneca-NGS/VarDictJava

if 'vardict' or 'vardict-java' is specified in the configuration. To use the
VarDict perl version:

https://github.com/AstraZeneca-NGS/VarDict

specify 'vardict-perl'.
"""
import os
import itertools
import sys

import toolz as tz
import pybedtools

from bcbio import bam, broad, utils
from bcbio.bam import highdepth
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils, shared
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import annotation, bamprep, vcfutils

def _is_bed_file(target):
    return target and isinstance(target, basestring) and os.path.isfile(target)

def _vardict_options_from_config(items, config, out_file, target=None):
    opts = ["-c 1", "-S 2", "-E 3", "-g 4"]
    # ["-z", "-F", "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-x", "0",
    #  "-k", "3", "-r", "4", "-m", "8"]

    resources = config_utils.get_resources("vardict", config)
    if resources.get("options"):
        opts += resources["options"]
    assert _is_bed_file(target)
    if any(tz.get_in(["config", "algorithm", "coverage_interval"], x, "").lower() == "genome"
            for x in items):
        target = shared.remove_highdepth_regions(target, items)
        target = shared.remove_lcr_regions(target, items)
    target = _enforce_max_region_size(target, items[0])
    opts += [target]  # this must be the last option
    return opts

def _enforce_max_region_size(in_file, data):
    """Ensure we don't have any chunks in the region greater than 1Mb.

    Larger sections have high memory usage on VarDictJava and failures
    on VarDict. This creates minimum windows from the input BED file
    to avoid these issues. Downstream VarDict merging sorts out any
    variants across windows.
    """
    max_size = 1e6
    overlap_size = 250
    def _has_larger_regions(f):
        return any(r.stop - r.start > max_size for r in pybedtools.BedTool(f))
    out_file = "%s-regionlimit%s" % utils.splitext_plus(in_file)
    if not utils.file_exists(out_file):
        if _has_larger_regions(in_file):
            with file_transaction(data, out_file) as tx_out_file:
                pybedtools.BedTool().window_maker(w=max_size,
                                                  s=max_size - overlap_size,
                                                  b=pybedtools.BedTool(in_file)).saveas(tx_out_file)
        else:
            utils.symlink_plus(in_file, out_file)
    return out_file

def run_vardict(align_bams, items, ref_file, assoc_files, region=None,
                  out_file=None):
    """Run VarDict variant calling.
    """
    if vcfutils.is_paired_analysis(align_bams, items):
        call_file = _run_vardict_paired(align_bams, items, ref_file,
                                        assoc_files, region, out_file)
    else:
        vcfutils.check_paired_problems(items)
        call_file = _run_vardict_caller(align_bams, items, ref_file,
                                        assoc_files, region, out_file)
    return call_file

def _get_jvm_opts(data, out_file):
    """Retrieve JVM options when running the Java version of VarDict.
    """
    if not dd.get_variantcaller(data).endswith("-perl"):
        resources = config_utils.get_resources("vardict", data["config"])
        jvm_opts = resources.get("jvm_opts", ["-Xms750m", "-Xmx4g"])
        jvm_opts += broad.get_default_jvm_opts(os.path.dirname(out_file))
        return "export VAR_DICT_OPTS='%s' && " % " ".join(jvm_opts)
    else:
        return ""

def _run_vardict_caller(align_bams, items, ref_file, assoc_files,
                          region=None, out_file=None):
    """Detect SNPs and indels with VarDict.
    """
    config = items[0]["config"]
    if out_file is None:
        out_file = "%s-variants.vcf.gz" % os.path.splitext(align_bams[0])[0]
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            target = shared.subset_variant_regions(dd.get_variant_regions(items[0]), region,
                                                   out_file, do_merge=False)
            for align_bam in align_bams:
                bam.index(align_bam, config)
            num_bams = len(align_bams)
            sample_vcf_names = []  # for individual sample names, given batch calling may be required
            for bamfile, item in itertools.izip(align_bams, items):
                # prepare commands
                sample = dd.get_sample_name(item)
                vardict = dd.get_variantcaller(items[0])
                vardict = "vardict-java" if not vardict.endswith("-perl") else "vardict"
                strandbias = "teststrandbias.R"
                var2vcf = "var2vcf_valid.pl"
                opts = (" ".join(_vardict_options_from_config(items, config, out_file, target))
                        if _is_bed_file(target) else "")
                vcfstreamsort = config_utils.get_program("vcfstreamsort", config)
                compress_cmd = "| bgzip -c" if out_file.endswith("gz") else ""
                freq = float(utils.get_in(config, ("algorithm", "min_allele_fraction"), 10)) / 100.0
                coverage_interval = utils.get_in(config, ("algorithm", "coverage_interval"), "exome")
                # for deep targeted panels, require 50 worth of coverage
                var2vcf_opts = " -v 50 " if highdepth.get_median_coverage(items[0]) > 5000 else ""
                fix_ambig = vcfutils.fix_ambiguous_cl()
                remove_dup = vcfutils.remove_dup_cl()
                jvm_opts = _get_jvm_opts(items[0], tx_out_file)
                cmd = ("{jvm_opts}{vardict} -G {ref_file} -f {freq} "
                        "-N {sample} -b {bamfile} {opts} "
                        "| {strandbias}"
                        "| {var2vcf} -N {sample} -E -f {freq} {var2vcf_opts} "
                        "| {fix_ambig} | {remove_dup} | {vcfstreamsort} {compress_cmd}")
                if num_bams > 1:
                    temp_file_prefix = out_file.replace(".gz", "").replace(".vcf", "") + item["name"][1]
                    tmp_out = temp_file_prefix + ".temp.vcf"
                    tmp_out += ".gz" if out_file.endswith("gz") else ""
                    sample_vcf_names.append(tmp_out)
                    with file_transaction(item, tmp_out) as tx_tmp_file:
                        if not _is_bed_file(target):
                            vcfutils.write_empty_vcf(tx_tmp_file, config, samples=[sample])
                        else:
                            cmd += " > {tx_tmp_file}"
                            do.run(cmd.format(**locals()), "Genotyping with VarDict: Inference", {})
                else:
                    if not _is_bed_file(target):
                        vcfutils.write_empty_vcf(tx_out_file, config, samples=[sample])
                    else:
                        cmd += " > {tx_out_file}"
                        do.run(cmd.format(**locals()), "Genotyping with VarDict: Inference", {})
            if num_bams > 1:
                # N.B. merge_variant_files wants region in 1-based end-inclusive
                # coordinates. Thus use bamprep.region_to_gatk
                vcfutils.merge_variant_files(orig_files=sample_vcf_names,
                                                out_file=tx_out_file, ref_file=ref_file,
                                                config=config, region=bamprep.region_to_gatk(region))
    out_file = (annotation.add_dbsnp(out_file, assoc_files["dbsnp"], config)
                if assoc_files.get("dbsnp") else out_file)
    return out_file

def _run_vardict_paired(align_bams, items, ref_file, assoc_files,
                          region=None, out_file=None):
    """Detect variants with Vardict.

    This is used for paired tumor / normal samples.
    """
    config = items[0]["config"]
    if out_file is None:
        out_file = "%s-paired-variants.vcf.gz" % os.path.splitext(align_bams[0])[0]
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            target = shared.subset_variant_regions(dd.get_variant_regions(items[0]), region,
                                                   out_file, do_merge=True)
            paired = vcfutils.get_paired_bams(align_bams, items)
            if not _is_bed_file(target):
                vcfutils.write_empty_vcf(tx_out_file, config,
                                         samples=[x for x in [paired.tumor_name, paired.normal_name] if x])
            else:
                if not paired.normal_bam:
                    ann_file = _run_vardict_caller(align_bams, items, ref_file,
                                                   assoc_files, region, out_file)
                    return ann_file
                vcffilter = config_utils.get_program("vcffilter", config)
                vardict = dd.get_variantcaller(items[0])
                vardict = "vardict-java" if not vardict.endswith("-perl") else "vardict"
                vcfstreamsort = config_utils.get_program("vcfstreamsort", config)
                strandbias = "testsomatic.R"
                var2vcf = "var2vcf_paired.pl"
                compress_cmd = "| bgzip -c" if out_file.endswith("gz") else ""
                freq = float(utils.get_in(config, ("algorithm", "min_allele_fraction"), 10)) / 100.0
                # merge bed file regions as amplicon VarDict is only supported in single sample mode
                opts = " ".join(_vardict_options_from_config(items, config, out_file, target))
                coverage_interval = utils.get_in(config, ("algorithm", "coverage_interval"), "exome")
                # for deep targeted panels, require 50 worth of coverage
                var2vcf_opts = " -v 50 " if highdepth.get_median_coverage(items[0]) > 5000 else ""
                fix_ambig = vcfutils.fix_ambiguous_cl()
                remove_dup = vcfutils.remove_dup_cl()
                if any("vardict_somatic_filter" in tz.get_in(("config", "algorithm", "tools_off"), data, [])
                       for data in items):
                    somatic_filter = ""
                else:
                    somatic_filter = ("| %s -x 'bcbio.variation.freebayes.call_somatic(x)'" %
                                      os.path.join(os.path.dirname(sys.executable), "py"))
                jvm_opts = _get_jvm_opts(items[0], tx_out_file)
                cmd = ("{jvm_opts}{vardict} -G {ref_file} -f {freq} "
                       "-N {paired.tumor_name} -b \"{paired.tumor_bam}|{paired.normal_bam}\" {opts} "
                       "| {strandbias} "
                       "| {var2vcf} -N \"{paired.tumor_name}|{paired.normal_name}\" -f {freq} {var2vcf_opts} "
                       "| bcftools filter -m '+' -s 'REJECT' -e 'STATUS !~ \".*Somatic\"' 2> /dev/null "
                       "| sed 's/\\\\.*Somatic\\\\/Somatic/' "
                       "| sed 's/REJECT,Description=\".*\">/REJECT,Description=\"Not Somatic via VarDict\">/' "
                       "{somatic_filter} | {fix_ambig} | {remove_dup} | {vcfstreamsort} "
                       "{compress_cmd} > {tx_out_file}")
                bam.index(paired.tumor_bam, config)
                bam.index(paired.normal_bam, config)
                do.run(cmd.format(**locals()), "Genotyping with VarDict: Inference", {})
    out_file = (annotation.add_dbsnp(out_file, assoc_files["dbsnp"], config)
                if assoc_files.get("dbsnp") else out_file)
    return out_file
