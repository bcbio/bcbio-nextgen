"""Bayesian variant calling with FreeBayes.

https://github.com/ekg/freebayes
"""

from distutils.version import LooseVersion
import os
import subprocess
import sys

import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.heterogeneity import chromhacks
from bcbio.pipeline import config_utils, shared
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import annotation, bedutils, ploidy, vcfutils
from bcbio.variation.vcfutils import (get_paired_bams, is_paired_analysis,
                                      move_vcf)

def region_to_freebayes(region):
    if isinstance(region, (list, tuple)):
        chrom, start, end = region
        return "%s:%s..%s" % (chrom, start, end)
    else:
        return region

def _freebayes_options_from_config(items, config, out_file, region=None):
    """Prepare standard options from configuration input.

    Input BED target files are merged to avoid overlapping regions which
    cause FreeBayes to call multiple times.

    Checks for empty sets of target regions after filtering for high depth,
    in which case we should skip the FreeBayes run.
    """
    opts = ["--genotype-qualities", "--strict-vcf"]
    cur_ploidy = ploidy.get_ploidy(items, region)
    base_ploidy = ploidy.get_ploidy(items)
    opts += ["--ploidy", str(cur_ploidy)]
    # Adjust min fraction when trying to call more sensitively in certain
    # regions. This is primarily meant for pooled mitochondrial calling.
    if (isinstance(region, (list, tuple)) and chromhacks.is_mitochondrial(region[0])
          and cur_ploidy >= base_ploidy and "--min-alternate-fraction" not in opts and "-F" not in opts):
        opts += ["--min-alternate-fraction", "0.01"]
    variant_regions = bedutils.merge_overlaps(bedutils.population_variant_regions(items), items[0])
    # Produce gVCF output
    if any("gvcf" in dd.get_tools_on(d) for d in items):
        opts += ["--gvcf", "--gvcf-chunk", "50000"]
    no_target_regions = False
    target = shared.subset_variant_regions(variant_regions, region, out_file, items)
    if target:
        if isinstance(target, basestring) and os.path.isfile(target):
            if any(tz.get_in(["config", "algorithm", "coverage_interval"], x, "").lower() == "genome"
                   for x in items):
                target = shared.remove_highdepth_regions(target, items)
                if os.path.getsize(target) == 0:
                    no_target_regions = True
            opts += ["--targets", target]
        else:
            opts += ["--region", region_to_freebayes(target)]
    resources = config_utils.get_resources("freebayes", config)
    if resources.get("options"):
        opts += resources["options"]
    return opts, no_target_regions

def _add_somatic_opts(opts, paired):
    """Add somatic options to current set. See _run_freebayes_paired for references.
    """
    if "--min-alternate-fraction" not in opts and "-F" not in opts:
        # add minimum reportable allele frequency
        # FreeBayes defaults to 20%, but use 10% by default for the
        # tumor case
        min_af = float(utils.get_in(paired.tumor_config, ("algorithm",
                                                          "min_allele_fraction"), 10)) / 100.0
        opts += " --min-alternate-fraction %s" % min_af
    # Recommended settings for cancer calling
    opts += (" --pooled-discrete --pooled-continuous "
             "--report-genotype-likelihood-max --allele-balance-priors-off")
    return opts

def run_freebayes(align_bams, items, ref_file, assoc_files, region=None,
                  out_file=None):
    """Run FreeBayes variant calling, either paired tumor/normal or germline calling.
    """
    if is_paired_analysis(align_bams, items):
        paired = get_paired_bams(align_bams, items)
        if not paired.normal_bam:
            call_file = _run_freebayes_caller(align_bams, items, ref_file,
                                              assoc_files, region, out_file, somatic=paired)
        else:
            call_file = _run_freebayes_paired(align_bams, items, ref_file,
                                              assoc_files, region, out_file)
    else:
        vcfutils.check_paired_problems(items)
        call_file = _run_freebayes_caller(align_bams, items, ref_file,
                                          assoc_files, region, out_file)

    return call_file

def _clean_freebayes_fmt_cl():
    """For FreeBayes pre 1.1.0, need to remove problematic DPR FORMAT annotation.

    This can be removed after bcbio 1.0.1 release which will default to
    FreeBayes 1.1.0
    """
    version = subprocess.check_output(["freebayes", "--version"])
    version = LooseVersion(version.split("version:")[-1].strip().replace("v", ""))
    return "bcftools annotate -x FMT/DPR |" if version < LooseVersion("1.1.0") else ""

def _run_freebayes_caller(align_bams, items, ref_file, assoc_files,
                          region=None, out_file=None, somatic=None):
    """Detect SNPs and indels with FreeBayes.

    Performs post-filtering to remove very low quality variants which
    can cause issues feeding into GATK. Breaks variants into individual
    allelic primitives for analysis and evaluation.
    """
    config = items[0]["config"]
    if out_file is None:
        out_file = "%s-variants.vcf.gz" % os.path.splitext(align_bams[0])[0]
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            freebayes = config_utils.get_program("freebayes", config)
            input_bams = " ".join("-b %s" % x for x in align_bams)
            opts, no_target_regions = _freebayes_options_from_config(items, config, out_file, region)
            if no_target_regions:
                vcfutils.write_empty_vcf(tx_out_file, config, samples=[dd.get_sample_name(d) for d in items])
            else:
                opts = " ".join(opts)
                # Recommended options from 1000 genomes low-complexity evaluation
                # https://groups.google.com/d/msg/freebayes/GvxIzjcpbas/1G6e3ArxQ4cJ
                opts += " --min-repeat-entropy 1"
                # Remove partial observations, which cause a preference for heterozygote calls
                # https://github.com/ekg/freebayes/issues/234#issuecomment-205331765
                opts += " --no-partial-observations"
                if somatic:
                    opts = _add_somatic_opts(opts, somatic)
                compress_cmd = "| bgzip -c" if out_file.endswith("gz") else ""
                fix_ambig = vcfutils.fix_ambiguous_cl()
                clean_fmt_cmd = _clean_freebayes_fmt_cl()
                py_cl = os.path.join(os.path.dirname(sys.executable), "py")
                cmd = ("{freebayes} -f {ref_file} {opts} {input_bams} "
                       """| bcftools filter -i 'ALT="<*>" || QUAL > 5' """
                       "| {fix_ambig} | {clean_fmt_cmd} bcftools view -a - | "
                       "{py_cl} -x 'bcbio.variation.freebayes.remove_missingalt(x)' | "
                       "vcfallelicprimitives -t DECOMPOSED --keep-geno | vcffixup - | vcfstreamsort | "
                       "vt normalize -n -r {ref_file} -q - | vcfuniqalleles | vt uniq - 2> /dev/null "
                       "{compress_cmd} > {tx_out_file}")
                do.run(cmd.format(**locals()), "Genotyping with FreeBayes", {})
    ann_file = annotation.annotate_nongatk_vcf(out_file, align_bams,
                                               assoc_files.get("dbsnp"),
                                               ref_file, config)
    return ann_file

def _run_freebayes_paired(align_bams, items, ref_file, assoc_files,
                          region=None, out_file=None):
    """Detect SNPs and indels with FreeBayes for paired tumor/normal samples.

    Sources of options for FreeBayes:
    mailing list: https://groups.google.com/d/msg/freebayes/dTWBtLyM4Vs/HAK_ZhJHguMJ
    mailing list: https://groups.google.com/forum/#!msg/freebayes/LLH7ZfZlVNs/63FdD31rrfEJ
    speedseq: https://github.com/cc2qe/speedseq/blob/e6729aa2589eca4e3a946f398c1a2bdc15a7300d/bin/speedseq#L916
    sga/freebayes: https://github.com/jts/sga-extra/blob/7e28caf71e8107b697f9be7162050e4fa259694b/
                   sga_generate_varcall_makefile.pl#L299
    """
    config = items[0]["config"]
    if out_file is None:
        out_file = "%s-paired-variants.vcf.gz" % os.path.splitext(align_bams[0])[0]
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            paired = get_paired_bams(align_bams, items)
            assert paired.normal_bam, "Require normal BAM for FreeBayes paired calling and filtering"

            freebayes = config_utils.get_program("freebayes", config)
            opts, no_target_regions = _freebayes_options_from_config(items, config, out_file, region)
            if no_target_regions:
                vcfutils.write_empty_vcf(tx_out_file, config,
                                         samples=[x for x in [paired.tumor_name, paired.normal_name] if x])
            else:
                opts = " ".join(opts)
                opts += " --min-repeat-entropy 1"
                opts += " --no-partial-observations"
                opts = _add_somatic_opts(opts, paired)
                compress_cmd = "| bgzip -c" if out_file.endswith("gz") else ""
                fix_ambig = vcfutils.fix_ambiguous_cl()
                clean_fmt_cmd = _clean_freebayes_fmt_cl()
                py_cl = os.path.join(os.path.dirname(sys.executable), "py")
                cl = ("{freebayes} -f {ref_file} {opts} "
                      "{paired.tumor_bam} {paired.normal_bam} "
                      """| bcftools filter -i 'ALT="<*>" || QUAL > 5' """
                      "| {py_cl} -x 'bcbio.variation.freebayes.call_somatic(x)' "
                      "| {fix_ambig} | {clean_fmt_cmd} bcftools view -a - | "
                      "{py_cl} -x 'bcbio.variation.freebayes.remove_missingalt(x)' | "
                      "vcfallelicprimitives -t DECOMPOSED --keep-geno | vcffixup - | vcfstreamsort | "
                      "vt normalize -n -r {ref_file} -q - | vcfuniqalleles | vt uniq - 2> /dev/null "
                      "{compress_cmd} > {tx_out_file}")
                do.run(cl.format(**locals()), "Genotyping paired variants with FreeBayes", {})
    ann_file = annotation.annotate_nongatk_vcf(out_file, align_bams,
                                               assoc_files.get("dbsnp"), ref_file,
                                               config)
    return ann_file

# ## Filtering

def _check_lods(parts, tumor_thresh, normal_thresh):
    """Ensure likelihoods for tumor and normal pass thresholds.

    Skipped if no FreeBayes GL annotations available.
    """
    try:
        gl_index = parts[8].split(":").index("GL")
    except ValueError:
        return True
    try:
        tumor_gls = [float(x) for x in parts[9].split(":")[gl_index].split(",") if x != "."]
        if tumor_gls:
            tumor_lod = max(tumor_gls[i] - tumor_gls[0] for i in range(1, len(tumor_gls)))
        else:
            tumor_lod = -1.0
    # No GL information, no tumor call (so fail it)
    except IndexError:
        tumor_lod = -1.0
    try:
        normal_gls = [float(x) for x in parts[10].split(":")[gl_index].split(",") if x != "."]
        if normal_gls:
            normal_lod = min(normal_gls[0] - normal_gls[i] for i in range(1, len(normal_gls)))
        else:
            normal_lod = normal_thresh
    # No GL inofmration, no normal call (so pass it)
    except IndexError:
        normal_lod = normal_thresh
    return normal_lod >= normal_thresh and tumor_lod >= tumor_thresh

def _check_freqs(parts):
    """Ensure frequency of tumor to normal passes a reasonable threshold.

    Avoids calling low frequency tumors also present at low frequency in normals,
    which indicates a contamination or persistent error.
    """
    thresh_ratio = 2.7
    try:  # FreeBayes
        ao_index = parts[8].split(":").index("AO")
        ro_index = parts[8].split(":").index("RO")
    except ValueError:
        ao_index, ro_index = None, None
    try:  # VarDict
        af_index = parts[8].split(":").index("AF")
    except ValueError:
        af_index = None
    if af_index is None and ao_index is None:
        # okay to skip if a gVCF record
        if parts[4].find("<*>") == -1:
            raise NotImplementedError("Unexpected format annotations: %s" % parts[8])
    def _calc_freq(item):
        try:
            if ao_index is not None and ro_index is not None:
                ao = sum([int(x) for x in item.split(":")[ao_index].split(",")])
                ro = int(item.split(":")[ro_index])
                freq = ao / float(ao + ro)
            elif af_index is not None:
                freq = float(item.split(":")[af_index])
            else:
                freq = 0.0
        except (IndexError, ValueError, ZeroDivisionError):
            freq = 0.0
        return freq
    tumor_freq, normal_freq = _calc_freq(parts[9]), _calc_freq(parts[10])
    return normal_freq <= 0.001 or normal_freq <= tumor_freq / thresh_ratio

def remove_missingalt(line):
    """Remove lines that are missing an alternative allele.

    During cleanup of extra alleles, bcftools has an issue in complicated cases
    with duplicate alleles and will end up stripping all alternative alleles.
    This removes those lines to avoid issues downstream.
    """
    if not line.startswith("#"):
        parts = line.split("\t")
        if parts[4] == ".":
            return None
    return line

def call_somatic(line):
    """Call SOMATIC variants from tumor/normal calls, adding REJECT filters and SOMATIC flag.

    Assumes tumor/normal called with tumor first and normal second, as done in bcbio
    implementation.

    Uses MuTect like somatic filter based on implementation in speedseq:
    https://github.com/cc2qe/speedseq/blob/e6729aa2589eca4e3a946f398c1a2bdc15a7300d/bin/speedseq#L62

    Extracts the genotype likelihoods (GLs) from FreeBayes, which are like phred scores
    except not multiplied by 10.0 (https://en.wikipedia.org/wiki/Phred_quality_score).
    For tumors, we retrieve the best likelihood to not be reference (the first GL) and
    for normal, the best likelhood to be reference.

    After calculating the likelihoods, we compare these to thresholds to pass variants
    at tuned sensitivity/precision. Tuning done on DREAM synthetic 3 dataset evaluations.

    We also check that the frequency of the tumor exceeds the frequency of the normal by
    a threshold to avoid calls that are low frequency in both tumor and normal. This supports
    both FreeBayes and VarDict output frequencies.
    """
    # Thresholds are like phred scores, so 3.5 = phred35
    tumor_thresh, normal_thresh = 3.5, 3.5
    if line.startswith("#CHROM"):
        headers = ['##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic event">',
                   ('##FILTER=<ID=REJECT,Description="Not somatic due to normal call frequency '
                    'or phred likelihoods: tumor: %s, normal %s.">')
                    % (int(tumor_thresh * 10), int(normal_thresh * 10))]
        return "\n".join(headers) + "\n" + line
    elif line.startswith("#"):
        return line
    else:
        parts = line.split("\t")
        if _check_lods(parts, tumor_thresh, normal_thresh) and _check_freqs(parts):
            parts[7] = parts[7] + ";SOMATIC"
        else:
            if parts[6] in set([".", "PASS"]):
                parts[6] = "REJECT"
            else:
                parts[6] += ";REJECT"
        line = "\t".join(parts)
        return line

def _clean_freebayes_output(line):
    """Clean FreeBayes output to make post-processing with GATK happy.

    XXX Not applied on recent versions which fix issues to be more compatible
    with bgzip output, but retained in case of need.

    - Remove lines from FreeBayes outputs where REF/ALT are identical:
      2       22816178        .       G       G       0.0339196
      or there are multiple duplicate alleles:
      4       60594753        .       TGAAA   T,T
    - Remove Type=Int specifications which are not valid VCF and GATK chokes
      on.
    """
    if line.startswith("#"):
        line = line.replace("Type=Int,D", "Type=Integer,D")
        return line
    else:
        parts = line.split("\t")
        alleles = [x.strip() for x in parts[4].split(",")] + [parts[3].strip()]
        if len(alleles) == len(set(alleles)):
            return line
    return None

def clean_vcf_output(orig_file, clean_fn, config, name="clean"):
    """Provide framework to clean a file in-place, with the specified clean
    function.
    """
    base, ext = utils.splitext_plus(orig_file)
    out_file = "{0}-{1}{2}".format(base, name, ext)
    if not utils.file_exists(out_file):
        with open(orig_file) as in_handle:
            with file_transaction(config, out_file) as tx_out_file:
                with open(tx_out_file, "w") as out_handle:
                    for line in in_handle:
                        update_line = clean_fn(line)
                        if update_line:
                            out_handle.write(update_line)
        move_vcf(orig_file, "{0}.orig".format(orig_file))
        move_vcf(out_file, orig_file)
        with open(out_file, "w") as out_handle:
            out_handle.write("Moved to {0}".format(orig_file))
