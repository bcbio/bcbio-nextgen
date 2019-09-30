"""Sensitive variant calling using VarDict.

Defaults to using the faster, equally sensitive Java port:

https://github.com/AstraZeneca-NGS/VarDictJava

if 'vardict' or 'vardict-java' is specified in the configuration. To use the
VarDict perl version:

https://github.com/AstraZeneca-NGS/VarDict

specify 'vardict-perl'.
"""
from decimal import *

from distutils.version import LooseVersion
import os
import sys
from six.moves import zip

import six
import toolz as tz
import pybedtools

from bcbio import bam, broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.heterogeneity import chromhacks
from bcbio.pipeline import config_utils, shared
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do, programs
from bcbio.variation import bamprep, bedutils, vcfutils

def _is_bed_file(target):
    return target and isinstance(target, six.string_types) and os.path.isfile(target)

def _vardict_options_from_config(items, config, out_file, target=None, is_rnaseq=False):
    var2vcf_opts = []
    opts = ["-c 1", "-S 2", "-E 3", "-g 4"]
    # ["-z", "-F", "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-x", "0",
    #  "-k", "3", "-r", "4", "-m", "8"]
    cores = dd.get_num_cores(items[0])
    if cores and cores > 1:
        opts += ["-th", str(cores)]
    # Disable SV calling for vardict, causes issues with regional analysis
    # by detecting SVs outside of target regions, which messes up merging
    # SV calling will be worked on as a separate step
    vardict_cl = get_vardict_command(items[0])
    version = programs.get_version_manifest(vardict_cl)
    if (vardict_cl and version and
        ((vardict_cl == "vardict-java" and LooseVersion(version) >= LooseVersion("1.5.5")) or
         (vardict_cl == "vardict" and LooseVersion(version) >= LooseVersion("2018.07.25")))):
        opts += ["--nosv"]
    if (vardict_cl and version and
         (vardict_cl == "vardict-java" and LooseVersion(version) >= LooseVersion("1.5.6"))):
        opts += ["--deldupvar"]
    # remove low mapping quality reads
    if not is_rnaseq:
        opts += ["-Q", "10"]
    # Remove QCfail reads, avoiding high depth repetitive regions
    opts += ["-F", "0x700"]
    resources = config_utils.get_resources("vardict", config)
    if resources.get("options"):
        opts += [str(x) for x in resources["options"]]
    resources = config_utils.get_resources("var2vcf", config)
    if resources.get("options"):
        var2vcf_opts += [str(x) for x in resources["options"]]
    if target and _is_bed_file(target):
        target = _enforce_max_region_size(target, items[0])
        opts += [target]  # this must be the last option
    _add_freq_options(config, opts, var2vcf_opts)
    return " ".join(opts), " ".join(var2vcf_opts)

def _add_freq_options(config, opts, var2vcf_opts):
    """ Setting -f option for vardict and var2vcf_valid
        Prioritizing settings in resources/vardict/options, then algorithm/min_allele_fraction:
    min_allele_fraction   "-f" in opts  var2vcfopts   ->   vardict -f            var2vcf -f
    yes                           yes   yes                opts                  var2vcfopts
    yes                           yes   -                  opts                  -
    yes                           -     yes                min_allele_fraction   var2vcfopts
    yes                           -     -                  min_allele_fraction   min_allele_fraction
    default                       yes   yes                opts                  var2vcfopts
    default                       yes   -                  opts                  -
    default                       -     yes                min_allele_fraction   var2vcfopts
    default                       -     -                  min_allele_fraction   min_allele_fraction
    """
    if "-f" not in opts:
        freq = Decimal(utils.get_in(config, ("algorithm", "min_allele_fraction"), 10)) / Decimal(100.0)
        opts.extend(["-f", str(freq)])
        if "-f" not in var2vcf_opts:
            var2vcf_opts.extend(["-f", str(freq)])

def _enforce_max_region_size(in_file, data):
    """Ensure we don't have any chunks in the region greater than 20kb.

    VarDict memory usage depends on size of individual windows in the input
    file. This breaks regions into 20kb chunks with 250bp overlaps. 20kb gives
    ~1Gb/core memory usage and the overlaps avoid missing indels spanning a
    gap. Downstream VarDict merging sorts out any variants across windows.

    https://github.com/AstraZeneca-NGS/VarDictJava/issues/64
    """
    max_size = 20000
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
    items = shared.add_highdepth_genome_exclusion(items)
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

    var2vcf_valid uses -A flag which reports all alleles and improves sensitivity:
    https://github.com/AstraZeneca-NGS/VarDict/issues/35#issuecomment-276738191
    """
    config = items[0]["config"]
    if out_file is None:
        out_file = "%s-variants.vcf.gz" % os.path.splitext(align_bams[0])[0]
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            vrs = bedutils.population_variant_regions(items)
            target = shared.subset_variant_regions(
                vrs, region, out_file, items=items, do_merge=False)
            num_bams = len(align_bams)
            sample_vcf_names = []  # for individual sample names, given batch calling may be required
            for bamfile, item in zip(align_bams, items):
                # prepare commands
                sample = dd.get_sample_name(item)
                vardict = get_vardict_command(items[0])
                opts, var2vcf_opts = _vardict_options_from_config(items, config, out_file, target)
                vcfstreamsort = config_utils.get_program("vcfstreamsort", config)
                compress_cmd = "| bgzip -c" if tx_out_file.endswith("gz") else ""
                fix_ambig_ref = vcfutils.fix_ambiguous_cl()
                fix_ambig_alt = vcfutils.fix_ambiguous_cl(5)
                remove_dup = vcfutils.remove_dup_cl()
                py_cl = os.path.join(utils.get_bcbio_bin(), "py")
                jvm_opts = _get_jvm_opts(items[0], tx_out_file)
                setup = ("%s && unset JAVA_HOME &&" % utils.get_R_exports())
                contig_cl = vcfutils.add_contig_to_header_cl(ref_file, tx_out_file)
                lowfreq_filter = _lowfreq_linear_filter(0, False)
                cmd = ("{setup}{jvm_opts}{vardict} -G {ref_file} "
                       "-N {sample} -b {bamfile} {opts} "
                       "| teststrandbias.R "
                       "| var2vcf_valid.pl -A -N {sample} -E {var2vcf_opts} "
                       "| {contig_cl} | bcftools filter -i 'QUAL >= 0' | {lowfreq_filter} "
                       "| {fix_ambig_ref} | {fix_ambig_alt} | {remove_dup} | {vcfstreamsort} {compress_cmd}")
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
    return out_file

def _lowfreq_linear_filter(tumor_index, is_paired):
    """Linear classifier for removing low frequency false positives.

    Uses a logistic classifier based on 0.5% tumor only variants from the smcounter2 paper:

    https://github.com/bcbio/bcbio_validations/tree/master/somatic-lowfreq

    The classifier uses strand bias (SBF) and read mismatches (NM) and
    applies only for low frequency (<2%) and low depth (<30) variants.
    """
    if is_paired:
        sbf = "FORMAT/SBF[%s]" % tumor_index
        nm = "FORMAT/NM[%s]" % tumor_index
    else:
        sbf = "INFO/SBF"
        nm = "INFO/NM"
    cmd = ("""bcftools filter --soft-filter 'LowFreqBias' --mode '+' """
           """-e  'FORMAT/AF[{tumor_index}:0] < 0.02 && FORMAT/VD[{tumor_index}] < 30 """
           """&& {sbf} < 0.1 && {nm} >= 2.0'""")
    return cmd.format(**locals())

def add_db_germline_flag(line):
    """Adds a DB flag for Germline filters, allowing downstream compatibility with PureCN.
    """
    if line.startswith("#CHROM"):
        headers = ['##INFO=<ID=DB,Number=0,Type=Flag,Description="Likely germline variant">']
        return "\n".join(headers) + "\n" + line
    elif line.startswith("#"):
        return line
    else:
        parts = line.split("\t")
        if parts[7].find("STATUS=Germline") >= 0:
            parts[7] += ";DB"
        return "\t".join(parts)

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
        qual = utils.safe_to_float(parts[5])
        dp = utils.safe_to_float(sample_ft.get("DP"))
        af = utils.safe_to_float(sample_ft.get("AF"))
        nm = utils.safe_to_float(sample_ft.get("NM"))
        mq = utils.safe_to_float(sample_ft.get("MQ"))
        ssfs = [x for x in parts[7].split(";") if x.startswith("SSF=")]
        pval = utils.safe_to_float(ssfs[0].split("=")[-1] if ssfs else None)
        fname = None
        if not chromhacks.is_sex(parts[0]) and dp is not None and af is not None:
            if dp * af < 6:
                if aligner == "bwa" and nm is not None and mq is not None:
                    if (mq < 55.0 and nm > 1.0) or (mq < 60.0 and nm > 2.0):
                        fname = "LowAlleleDepth"
                if dp < 10:
                    fname = "LowAlleleDepth"
                if qual is not None and qual < 45:
                    fname = "LowAlleleDepth"
        if af is not None and qual is not None and pval is not None:
            if af < 0.2 and qual < 45 and pval > 0.06:
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
            vrs = bedutils.population_variant_regions(items)
            target = shared.subset_variant_regions(vrs, region,
                                                   out_file, items=items, do_merge=True)
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
                compress_cmd = "| bgzip -c" if out_file.endswith("gz") else ""
                freq = float(utils.get_in(config, ("algorithm", "min_allele_fraction"), 10)) / 100.0
                # merge bed file regions as amplicon VarDict is only supported in single sample mode
                opts, var2vcf_opts = _vardict_options_from_config(items, config, out_file, target)
                fix_ambig_ref = vcfutils.fix_ambiguous_cl()
                fix_ambig_alt = vcfutils.fix_ambiguous_cl(5)
                remove_dup = vcfutils.remove_dup_cl()
                if any("vardict_somatic_filter" in tz.get_in(("config", "algorithm", "tools_off"), data, [])
                       for data in items):
                    somatic_filter = ""
                    freq_filter = ""
                else:
                    var2vcf_opts += " -M "  # this makes VarDict soft filter non-differential variants
                    somatic_filter = ("| sed 's/\\\\.*Somatic\\\\/Somatic/' "
                                      "| sed 's/REJECT,Description=\".*\">/REJECT,Description=\"Not Somatic via VarDict\">/' "
                                      """| %s -c 'from bcbio.variation import freebayes; """
                                      """freebayes.call_somatic("%s", "%s")' """
                                      % (sys.executable, paired.tumor_name, paired.normal_name))
                    freq_filter = ("| bcftools filter -m '+' -s 'REJECT' -e 'STATUS !~ \".*Somatic\"' 2> /dev/null "
                                   "| %s -x 'bcbio.variation.vardict.add_db_germline_flag(x)' "
                                   "| %s "
                                   "| %s -x 'bcbio.variation.vardict.depth_freq_filter(x, %s, \"%s\")'" %
                                   (os.path.join(os.path.dirname(sys.executable), "py"),
                                    _lowfreq_linear_filter(0, True),
                                    os.path.join(os.path.dirname(sys.executable), "py"),
                                    0, bam.aligner_from_header(paired.tumor_bam)))
                jvm_opts = _get_jvm_opts(items[0], tx_out_file)
                py_cl = os.path.join(utils.get_bcbio_bin(), "py")
                setup = ("%s && unset JAVA_HOME &&" % utils.get_R_exports())
                contig_cl = vcfutils.add_contig_to_header_cl(ref_file, tx_out_file)
                cmd = ("{setup}{jvm_opts}{vardict} -G {ref_file} "
                       "-N {paired.tumor_name} -b \"{paired.tumor_bam}|{paired.normal_bam}\" {opts} "
                       "| awk 'NF>=48' | testsomatic.R "
                       "| var2vcf_paired.pl -P 0.9 -m 4.25 {var2vcf_opts} "
                       "-N \"{paired.tumor_name}|{paired.normal_name}\" "
                       "| {contig_cl} {freq_filter} "
                       "| bcftools filter -i 'QUAL >= 0' "
                       "{somatic_filter} | {fix_ambig_ref} | {fix_ambig_alt} | {remove_dup} | {vcfstreamsort} "
                       "{compress_cmd} > {tx_out_file}")
                do.run(cmd.format(**locals()), "Genotyping with VarDict: Inference", {})
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
