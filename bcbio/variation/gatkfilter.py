"""Perform GATK based filtering: soft filters, CNNs and VQSR.
"""
import os
import gzip
from distutils.version import LooseVersion

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import vcfutils, vfilter

def run(call_file, ref_file, vrn_files, data):
    """Run filtering on the input call file, handling SNPs and indels separately.
    """
    algs = [data["config"]["algorithm"]] * len(data.get("vrn_files", [1]))
    if includes_missingalt(data):
        logger.info("Removing variants with missing alts from %s." % call_file)
        call_file = gatk_remove_missingalt(call_file, data)

    if "gatkcnn" in dd.get_tools_on(data):
        return _cnn_filter(call_file, vrn_files, data)
    elif config_utils.use_vqsr(algs, call_file):
        if vcfutils.is_gvcf_file(call_file):
            raise ValueError("Cannot force gVCF output with joint calling using tools_on: [gvcf] and use VQSR. "
                             "Try using cutoff-based soft filtering with tools_off: [vqsr]")
        snp_file, indel_file = vcfutils.split_snps_indels(call_file, ref_file, data["config"])
        snp_filter_file = _variant_filtration(snp_file, ref_file, vrn_files, data, "SNP",
                                              vfilter.gatk_snp_cutoff)
        indel_filter_file = _variant_filtration(indel_file, ref_file, vrn_files, data, "INDEL",
                                                vfilter.gatk_indel_cutoff)
        orig_files = [snp_filter_file, indel_filter_file]
        out_file = "%scombined.vcf.gz" % os.path.commonprefix(orig_files)
        combined_file = vcfutils.combine_variant_files(orig_files, out_file, ref_file, data["config"])
        return combined_file
    else:
        snp_filter = vfilter.gatk_snp_cutoff(call_file, data)
        indel_filter = vfilter.gatk_indel_cutoff(snp_filter, data)
        return indel_filter

# ## Convolutional Neural Networks (CNN)

def _cnn_filter(in_file, vrn_files, data):
    """Perform CNN filtering on input VCF using pre-trained models.
    """
    #tensor_type = "reference"  # 1D, reference sequence
    tensor_type = "read_tensor"  # 2D, reads, flags, mapping quality
    score_file = _cnn_score_variants(in_file, tensor_type, data)
    return _cnn_tranch_filtering(score_file, vrn_files, tensor_type, data)

def _cnn_tranch_filtering(in_file, vrn_files, tensor_type, data):
    """Filter CNN scored VCFs in tranches using standard SNP and Indel truth sets.
    """
    out_file = "%s-filter.vcf.gz" % utils.splitext_plus(in_file)[0]
    if not utils.file_uptodate(out_file, in_file):
        runner = broad.runner_from_config(data["config"])
        gatk_type = runner.gatk_type()
        assert gatk_type == "gatk4", "CNN filtering requires GATK4"
        if "train_hapmap" not in vrn_files:
            raise ValueError("CNN filtering requires HapMap training inputs: %s" % vrn_files)
        with file_transaction(data, out_file) as tx_out_file:
            params = ["-T", "FilterVariantTranches", "--variant", in_file,
                      "--output", tx_out_file,
                      "--snp-truth-vcf", vrn_files["train_hapmap"],
                      "--indel-truth-vcf", vrn_files["train_indels"]]
            if tensor_type == "reference":
                params += ["--info-key", "CNN_1D", "--tranche", "99"]
            else:
                assert tensor_type == "read_tensor"
                params += ["--info-key", "CNN_2D", "--tranche", "99"]
            runner.run_gatk(params)
    return vcfutils.bgzip_and_index(out_file, data["config"])

def _cnn_score_variants(in_file, tensor_type, data):
    """Score variants with pre-trained CNN models.
    """
    out_file = "%s-cnnscore.vcf.gz" % utils.splitext_plus(in_file)[0]
    if not utils.file_uptodate(out_file, in_file):
        runner = broad.runner_from_config(data["config"])
        gatk_type = runner.gatk_type()
        assert gatk_type == "gatk4", "CNN filtering requires GATK4"
        with file_transaction(data, out_file) as tx_out_file:
            params = ["-T", "CNNScoreVariants", "--variant", in_file, "--reference", dd.get_ref_file(data),
                    "--output", tx_out_file, "--input", dd.get_align_bam(data)]
            params += ["--tensor-type", tensor_type]
            runner.run_gatk(params)
    return vcfutils.bgzip_and_index(out_file, data["config"])

# ## Variant Quality Score Recalibration (VQSR)

def _apply_vqsr(in_file, ref_file, recal_file, tranch_file,
                sensitivity_cutoff, filter_type, data):
    """Apply VQSR based on the specified tranche, returning a filtered VCF file.
    """
    base, ext = utils.splitext_plus(in_file)
    out_file = "{base}-{filter}filter{ext}".format(base=base, ext=ext,
                                                   filter=filter_type)
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            broad_runner = broad.runner_from_config(data["config"])
            gatk_type = broad_runner.gatk_type()
            if gatk_type == "gatk4":
                params = ["-T", "ApplyVQSR",
                          "--variant", in_file,
                          "--output", tx_out_file,
                          "--recal-file", recal_file,
                          "--tranches-file", tranch_file]
            else:
                params = ["-T", "ApplyRecalibration",
                          "--input", in_file,
                          "--out", tx_out_file,
                          "--recal_file", recal_file,
                          "--tranches_file", tranch_file]
            params += ["-R", ref_file,
                       "--mode", filter_type]
            resources = config_utils.get_resources("gatk_apply_recalibration", data["config"])
            opts = resources.get("options", [])
            if not opts:
                if gatk_type == "gatk4":
                    opts += ["--truth-sensitivity-filter-level", sensitivity_cutoff]
                else:
                    opts += ["--ts_filter_level", sensitivity_cutoff]
            params += opts
            broad_runner.run_gatk(params)
    return out_file

def _get_training_data(vrn_files):
    """Retrieve training data, returning an empty set of information if not available.
    """
    out = {"SNP": [], "INDEL": []}
    # SNPs
    for name, train_info in [("train_hapmap", "known=false,training=true,truth=true,prior=15.0"),
                             ("train_omni", "known=false,training=true,truth=true,prior=12.0"),
                             ("train_1000g", "known=false,training=true,truth=false,prior=10.0"),
                             ("dbsnp", "known=true,training=false,truth=false,prior=2.0")]:
        if name not in vrn_files:
            return {}
        else:
            out["SNP"].append((name.replace("train_", ""), train_info, vrn_files[name]))
    # Indels
    if "train_indels" in vrn_files:
        out["INDEL"].append(("mills", "known=true,training=true,truth=true,prior=12.0",
                             vrn_files["train_indels"]))
    else:
        return {}
    return out

def _have_training_data(vrn_files):
    return len(_get_training_data(vrn_files)) > 0

def _get_vqsr_training(filter_type, vrn_files, gatk_type):
    """Return parameters for VQSR training, handling SNPs and Indels.
    """
    params = []
    for name, train_info, fname in _get_training_data(vrn_files)[filter_type]:
        if gatk_type == "gatk4":
            params.extend(["--resource:%s,%s" % (name, train_info), fname])
            if filter_type == "INDEL":
                params.extend(["--max-gaussians", "4"])
        else:
            params.extend(["-resource:%s,VCF,%s" % (name, train_info), fname])
            if filter_type == "INDEL":
                params.extend(["--maxGaussians", "4"])
    return params

def _get_vqsr_annotations(filter_type, data):
    """Retrieve appropriate annotations to use for VQSR based on filter type.

    Issues reported with MQ and bwa-mem quality distribution, results in intermittent
    failures to use VQSR:
    http://gatkforums.broadinstitute.org/discussion/4425/variant-recalibration-failing
    http://gatkforums.broadinstitute.org/discussion/4248/variantrecalibrator-removing-all-snps-from-the-training-set
    """
    if filter_type == "SNP":
        # MQ, MQRankSum
        anns = ["QD", "FS", "ReadPosRankSum", "SOR"]
    else:
        assert filter_type == "INDEL"
        # MQRankSum
        anns = ["QD", "FS", "ReadPosRankSum", "SOR"]
    if dd.get_coverage_interval(data) == "genome":
        anns += ["DP"]
    return anns

def _run_vqsr(in_file, ref_file, vrn_files, sensitivity_cutoff, filter_type, data):
    """Run variant quality score recalibration.
    """
    cutoffs = ["100.0", "99.99", "99.98", "99.97", "99.96", "99.95", "99.94", "99.93", "99.92", "99.91",
               "99.9", "99.8", "99.7", "99.6", "99.5", "99.0", "98.0", "90.0"]
    if sensitivity_cutoff not in cutoffs:
        cutoffs.append(sensitivity_cutoff)
        cutoffs.sort()
    broad_runner = broad.runner_from_config(data["config"])
    gatk_type = broad_runner.gatk_type()
    base = utils.splitext_plus(in_file)[0]
    recal_file = ("%s-vqsrrecal.vcf.gz" % base) if gatk_type == "gatk4" else ("%s.recal" % base)
    tranches_file = "%s.tranches" % base
    plot_file = "%s-plots.R" % base
    if not utils.file_exists(recal_file):
        with file_transaction(data, recal_file, tranches_file, plot_file) as (tx_recal, tx_tranches, tx_plot_file):
            params = ["-T", "VariantRecalibrator",
                      "-R", ref_file,
                      "--mode", filter_type]
            if gatk_type == "gatk4":
                params += ["--variant", in_file, "--output", tx_recal,
                           "--tranches-file", tx_tranches, "--rscript-file", tx_plot_file]
            else:
                params += ["--input", in_file, "--recal_file", tx_recal,
                           "--tranches_file", tx_tranches, "--rscript_file", tx_plot_file]
            params += _get_vqsr_training(filter_type, vrn_files, gatk_type)
            resources = config_utils.get_resources("gatk_variant_recalibrator", data["config"])
            opts = resources.get("options", [])
            if not opts:
                for cutoff in cutoffs:
                    opts += ["-tranche", str(cutoff)]
                for a in _get_vqsr_annotations(filter_type, data):
                    opts += ["-an", a]
            params += opts
            cores = dd.get_cores(data)
            memscale = {"magnitude": 0.9 * cores, "direction": "increase"} if cores > 1 else None
            try:
                broad_runner.new_resources("gatk-vqsr")
                broad_runner.run_gatk(params, log_error=False, memscale=memscale, parallel_gc=True)
            except:  # Can fail to run if not enough values are present to train.
                return None, None
    if gatk_type == "gatk4":
        vcfutils.bgzip_and_index(recal_file, data["config"])
    return recal_file, tranches_file

# ## SNP and indel specific variant filtration

def _already_cutoff_filtered(in_file, filter_type):
    """Check if we have a pre-existing cutoff-based filter file from previous VQSR failure.
    """
    filter_file = "%s-filter%s.vcf.gz" % (utils.splitext_plus(in_file)[0], filter_type)
    return utils.file_exists(filter_file)

def _variant_filtration(in_file, ref_file, vrn_files, data, filter_type,
                        hard_filter_fn):
    """Filter SNP and indel variant calls using GATK best practice recommendations.

    Use cutoff-based soft filters if configuration indicates too little data or
    already finished a cutoff-based filtering step, otherwise try VQSR.
    """
    # Algorithms multiplied by number of input files to check for large enough sample sizes
    algs = [data["config"]["algorithm"]] * len(data.get("vrn_files", [1]))
    if (not config_utils.use_vqsr(algs, in_file) or
          _already_cutoff_filtered(in_file, filter_type)):
        logger.info("Skipping VQSR, using cutoff-based filers: we don't have whole genome input data")
        return hard_filter_fn(in_file, data)
    elif not _have_training_data(vrn_files):
        logger.info("Skipping VQSR, using cutoff-based filers: genome build does not have sufficient training data")
        return hard_filter_fn(in_file, data)
    else:
        sensitivities = {"INDEL": "98.0", "SNP": "99.97"}
        recal_file, tranches_file = _run_vqsr(in_file, ref_file, vrn_files,
                                              sensitivities[filter_type], filter_type, data)
        if recal_file is None:  # VQSR failed
            logger.info("VQSR failed due to lack of training data. Using cutoff-based soft filtering.")
            return hard_filter_fn(in_file, data)
        else:
            return _apply_vqsr(in_file, ref_file, recal_file, tranches_file,
                               sensitivities[filter_type], filter_type, data)

def includes_missingalt(data):
    """
    As of GATK 4.1.0.0, variants with missing alts are generated
    (see https://github.com/broadinstitute/gatk/issues/5650)
    """
    MISSINGALT_VERSION = LooseVersion("4.1.0.0")
    version = LooseVersion(broad.get_gatk_version(config=dd.get_config(data)))
    return version >= MISSINGALT_VERSION

def gatk_remove_missingalt(in_file, data):
    """
    GATK 4.1.0.0 outputs variants that have missing ALTs, which breaks downstream
    tools, this filters those out.
    """
    base = in_file.split('.vcf.gz')[0]
    out_file = "%s-nomissingalt%s" % (base, '.vcf.gz')
    if utils.file_exists(out_file):
        return out_file
    no_gzip_out = out_file.replace(".vcf.gz", ".vcf")
    with file_transaction(no_gzip_out) as tx_out_file:
        with utils.open_gzipsafe(in_file) as in_handle, open(tx_out_file, "w") as out_handle:
            for line in in_handle:
                line = remove_missingalt(line)
                if line:
                    out_handle.write(line)
    return vcfutils.bgzip_and_index(no_gzip_out, data["config"])

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
