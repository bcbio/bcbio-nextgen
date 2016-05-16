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

from bcbio import broad, utils
from bcbio.bam import highdepth
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils, shared
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import annotation, bamprep, bedutils, vcfutils

def _is_bed_file(target):
    return target and isinstance(target, basestring) and os.path.isfile(target)

def _vardict_options_from_config(items, config, out_file, target=None):
    opts = ["-c 1", "-S 2", "-E 3", "-g 4"]
    # ["-z", "-F", "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-x", "0",
    #  "-k", "3", "-r", "4", "-m", "8"]

    # remove low mapping quality reads
    opts += ["-Q", "10"]
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
    if get_vardict_command(data) == "vardict-java":
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
            vrs = bedutils.population_variant_regions(items)
            target = shared.subset_variant_regions(vrs, region, out_file, do_merge=False)
            num_bams = len(align_bams)
            sample_vcf_names = []  # for individual sample names, given batch calling may be required
            for bamfile, item in itertools.izip(align_bams, items):
                # prepare commands
                sample = dd.get_sample_name(item)
                vardict = get_vardict_command(items[0])
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
                r_setup = "unset R_HOME && export PATH=%s:$PATH && " % os.path.dirname(utils.Rscript_cmd())
                cmd = ("{r_setup}{jvm_opts}{vardict} -G {ref_file} -f {freq} "
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

def _safe_to_float(x):
    if x is None:
        return None
    else:
        try:
            return float(x)
        except ValueError:
            return None

def depth_freq_filter(line, tumor_index, aligner):
    """Command line to filter VarDict calls based on depth, frequency and quality.

    Looks at regions with low depth for allele frequency (AF * DP < 6, the equivalent
    of < 13bp for heterogygote calls, but generalized. Within these calls filters if a
    calls has:

    - Low mapping quality and multiple mismatches in a read (NM)
        For bwa only: MQ < 55.0 and NM > 1.0 or MQ < 60.0 and NM > 2.0
    - Low depth (DP < 10)
    - Low QUAL (QUAL < 45)

    Also filters in low allele frequency regions with poor quality, if all of these are
    true:
    - Allele frequency < 0.2
    - Quality < 55
    - P-value (SSF) > 0.06
    """
    if line.startswith("#CHROM"):
        headers = [('##FILTER=<ID=LowAlleleDepth,Description="Low depth per allele frequency '
                    'along with poor depth, quality, mapping quality and read mismatches.">'),
                   ('##FILTER=<ID=LowFreqQuality,Description="Low frequency read with '
                    'poor quality and p-value (SSF).">')]
        return "\n".join(headers) + "\n" + line
    elif line.startswith("#"):
        return line
    else:
        parts = line.split("\t")
        sample_ft = {a: v for (a, v) in zip(parts[8].split(":"), parts[9 + tumor_index].split(":"))}
        qual = _safe_to_float(parts[5])
        dp = _safe_to_float(sample_ft.get("DP"))
        af = _safe_to_float(sample_ft.get("AF"))
        nm = _safe_to_float(sample_ft.get("NM"))
        mq = _safe_to_float(sample_ft.get("MQ"))
        ssfs = [x for x in parts[7].split(";") if x.startswith("SSF=")]
        pval = _safe_to_float(ssfs[0].split("=")[-1] if ssfs else None)
        fname = None
        if dp is not None and af is not None:
            if dp * af < 6:
                if aligner == "bwa" and nm is not None and mq is not None:
                    if (mq < 55.0 and nm > 1.0) or (mq < 60.0 and nm > 2.0):
                        fname = "LowAlleleDepth"
                if dp < 10:
                    fname = "LowAlleleDepth"
                if qual is not None and qual < 45:
                    fname = "LowAlleleDepth"
        if af is not None and qual is not None and pval is not None:
            if af < 0.2 and qual < 55 and pval > 0.06:
                fname = "LowFreqQuality"
        if fname:
            if parts[6] in set([".", "PASS"]):
                parts[6] = fname
            else:
                parts[6] += ";%s" % fname
        line = "\t".join(parts)
        return line

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
                vardict = get_vardict_command(items[0])
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
                    freq_filter = ""
                else:
                    var2vcf_opts += " -M "  # this makes VarDict soft filter non-differential variants
                    somatic_filter = ("| sed 's/\\\\.*Somatic\\\\/Somatic/' "
                                      "| sed 's/REJECT,Description=\".*\">/REJECT,Description=\"Not Somatic via VarDict\">/' "
                                      "| %s -x 'bcbio.variation.freebayes.call_somatic(x)'" %
                                      os.path.join(os.path.dirname(sys.executable), "py"))
                    freq_filter = ("| bcftools filter -m '+' -s 'REJECT' -e 'STATUS !~ \".*Somatic\"' 2> /dev/null "
                                   "| %s -x 'bcbio.variation.vardict.depth_freq_filter(x, %s, \"%s\")'" %
                                   (os.path.join(os.path.dirname(sys.executable), "py"),
                                     0, dd.get_aligner(paired.tumor_data)))
                jvm_opts = _get_jvm_opts(items[0], tx_out_file)
                r_setup = "unset R_HOME && export PATH=%s:$PATH && " % os.path.dirname(utils.Rscript_cmd())
                cmd = ("{r_setup}{jvm_opts}{vardict} -G {ref_file} -f {freq} "
                       "-N {paired.tumor_name} -b \"{paired.tumor_bam}|{paired.normal_bam}\" {opts} "
                       "| {strandbias} "
                       "| {var2vcf} -P 0.9 -m 4.25 -f {freq} {var2vcf_opts} "
                       "-N \"{paired.tumor_name}|{paired.normal_name}\" "
                       "{freq_filter} "
                       "{somatic_filter} | {fix_ambig} | {remove_dup} | {vcfstreamsort} "
                       "{compress_cmd} > {tx_out_file}")
                do.run(cmd.format(**locals()), "Genotyping with VarDict: Inference", {})
    out_file = (annotation.add_dbsnp(out_file, assoc_files["dbsnp"], config)
                if assoc_files.get("dbsnp") else out_file)
    return out_file

def get_vardict_command(data):
    """
    convert variantcaller specification to proper vardict command, handling
    string or list specification
    """
    vcaller = dd.get_variantcaller(data)
    if isinstance(vcaller, list):
        vardict = [x for x in vcaller if "vardict" in x]
        if not vardict:
            return None
        vardict = vardict[0]
    elif not vcaller:
        return None
    else:
        vardict = vcaller
    vardict = "vardict-java" if not vardict.endswith("-perl") else "vardict"
    return vardict
