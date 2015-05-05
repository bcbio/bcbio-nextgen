"""Perform GATK based filtering, perferring variant quality score recalibration.

Performs hard filtering when VQSR fails on smaller sets of variant calls.
"""
import os

import toolz as tz

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
    snp_file, indel_file = vcfutils.split_snps_indels(call_file, ref_file, data["config"])
    snp_filter_file = _variant_filtration(snp_file, ref_file, vrn_files, data, "SNP",
                                          vfilter.gatk_snp_hard)
    indel_filter_file = _variant_filtration(indel_file, ref_file, vrn_files, data, "INDEL",
                                            vfilter.gatk_indel_hard)
    orig_files = [snp_filter_file, indel_filter_file]
    out_file = "%scombined.vcf.gz" % os.path.commonprefix(orig_files)
    combined_file = vcfutils.combine_variant_files(orig_files, out_file, ref_file, data["config"])
    return _filter_nonref(combined_file, data)

def _filter_nonref(in_file, data):
    """Remove NON_REF gVCF items from GATK VCF output; these occasionally sneak through in joint calling.
    """
    out_file = "%s-gatkclean%s" % utils.splitext_plus(in_file)
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = "gunzip -c {in_file} | grep -v NON_REF | bgzip -c > {tx_out_file}"
            do.run(cmd.format(**locals()), "Remove stray NON_REF gVCF information from VCF output", data)
        vcfutils.bgzip_and_index(out_file, data["config"])
    return out_file

def _apply_vqsr(in_file, ref_file, recal_file, tranch_file,
                sensitivity_cutoff, filter_type, data):
    """Apply VQSR based on the specified tranche, returning a filtered VCF file.
    """
    broad_runner = broad.runner_from_config(data["config"])
    base, ext = utils.splitext_plus(in_file)
    out_file = "{base}-{filter}filter{ext}".format(base=base, ext=ext,
                                                   filter=filter_type)
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            params = ["-T", "ApplyRecalibration",
                      "-R", ref_file,
                      "--input", in_file,
                      "--out", tx_out_file,
                      "--ts_filter_level", sensitivity_cutoff,
                      "--tranches_file", tranch_file,
                      "--recal_file", recal_file,
                      "--mode", filter_type]
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

def _get_vqsr_training(filter_type, vrn_files):
    """Return parameters for VQSR training, handling SNPs and Indels.
    """
    params = []
    for name, train_info, fname in _get_training_data(vrn_files)[filter_type]:
        params.extend(["-resource:%s,VCF,%s" % (name, train_info), fname])
    if filter_type == "INDEL":
        params.extend(["--maxGaussians", "4"])
    return params

def _get_vqsr_annotations(filter_type):
    """Retrieve appropriate annotations to use for VQSR based on filter type.

    Issues reported with MQ and bwa-mem quality distribution, results in intermittent
    failures to use VQSR:
    http://gatkforums.broadinstitute.org/discussion/4425/variant-recalibration-failing
    http://gatkforums.broadinstitute.org/discussion/4248/variantrecalibrator-removing-all-snps-from-the-training-set
    """
    if filter_type == "SNP":
        # MQ, MQRankSum
        return ["DP", "QD", "FS", "ReadPosRankSum"]
    else:
        assert filter_type == "INDEL"
        # MQRankSum
        return ["DP", "QD", "FS", "ReadPosRankSum"]

def _run_vqsr(in_file, ref_file, vrn_files, sensitivity_cutoff, filter_type, data):
    """Run variant quality score recalibration.
    """
    cutoffs = ["100.0", "99.99", "99.98", "99.97", "99.96", "99.95", "99.94", "99.93", "99.92", "99.91",
               "99.9", "99.8", "99.7", "99.6", "99.5", "99.0", "98.0", "90.0"]
    if sensitivity_cutoff not in cutoffs:
        cutoffs.append(sensitivity_cutoff)
        cutoffs.sort()
    broad_runner = broad.runner_from_config(data["config"])
    base = utils.splitext_plus(in_file)[0]
    recal_file = "%s.recal" % base
    tranches_file = "%s.tranches" % base
    if not utils.file_exists(recal_file):
        with file_transaction(data, recal_file, tranches_file) as (tx_recal, tx_tranches):
            params = ["-T", "VariantRecalibrator",
                      "-R", ref_file,
                      "--input", in_file,
                      "--mode", filter_type,
                      "--recal_file", tx_recal,
                      "--tranches_file", tx_tranches]
            for cutoff in cutoffs:
                params += ["-tranche", str(cutoff)]
            params += _get_vqsr_training(filter_type, vrn_files)
            for a in _get_vqsr_annotations(filter_type):
                params += ["-an", a]
            cores = dd.get_cores(data)
            memscale = {"magnitude": 0.9 * cores, "direction": "increase"} if cores > 1 else None
            try:
                broad_runner.new_resources("gatk-vqsr")
                broad_runner.run_gatk(params, log_error=False, memscale=memscale)
            except:  # Can fail to run if not enough values are present to train.
                return None, None
    return recal_file, tranches_file

# ## SNP and indel specific variant filtration

def _already_hard_filtered(in_file, filter_type):
    """Check if we have a pre-existing hard filter file from previous VQSR failure.
    """
    filter_file = "%s-filter%s.vcf.gz" % (utils.splitext_plus(in_file)[0], filter_type)
    return utils.file_exists(filter_file)

def _variant_filtration(in_file, ref_file, vrn_files, data, filter_type,
                        hard_filter_fn):
    """Filter SNP and indel variant calls using GATK best practice recommendations.

    Hard filter if configuration indicates too little data or already finished a
    hard filtering, otherwise try VQSR.
    """
    # Algorithms multiplied by number of input files to check for large enough sample sizes
    algs = [data["config"]["algorithm"]] * len(data.get("vrn_files", [1]))
    if (not config_utils.use_vqsr(algs) or
          _already_hard_filtered(in_file, filter_type)):
        logger.info("Skipping VQSR, using hard filers: we don't have whole genome input data")
        return hard_filter_fn(in_file, data)
    elif not _have_training_data(vrn_files):
        logger.info("Skipping VQSR, using hard filers: genome build does not have sufficient training data")
        return hard_filter_fn(in_file, data)
    else:
        sensitivities = {"INDEL": "98.0", "SNP": "99.97"}
        recal_file, tranches_file = _run_vqsr(in_file, ref_file, vrn_files,
                                              sensitivities[filter_type], filter_type, data)
        if recal_file is None:  # VQSR failed
            logger.info("VQSR failed due to lack of training data. Using hard filtering.")
            return hard_filter_fn(in_file, data)
        else:
            return _apply_vqsr(in_file, ref_file, recal_file, tranches_file,
                               sensitivities[filter_type], filter_type, data)
