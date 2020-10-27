"""Ensemble methods that create consensus calls from multiple approaches.

This handles merging calls produced by multiple calling methods or
technologies into a single consolidated callset. Uses the bcbio.variation
toolkit: https://github.com/chapmanb/bcbio.variation and bcbio.variation.recall:
https://github.com/chapmanb/bcbio.variation.recall
"""
import collections
import copy
import glob
import math
import os

import yaml
import toolz as tz

from bcbio import utils
from bcbio.cwl import cwlutils
from bcbio.log import logger
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import population, validate, vcfutils, multi, normalize, effects

def batch(samples):
    """CWL: batch together per sample, joint and germline calls for ensemble combination.

    Sets up groups of same sample/batch variant calls for ensemble calling, as
    long as we have more than one caller per group.
    """
    samples = [utils.to_single_data(x) for x in samples]
    sample_order = [dd.get_sample_name(x) for x in samples]
    batch_groups = collections.defaultdict(list)
    for data in samples:
        batch_samples = tuple(data.get("batch_samples", [dd.get_sample_name(data)]))
        batch_groups[(batch_samples, dd.get_phenotype(data))].append(data)

    out = []
    for (batch_samples, phenotype), gsamples in batch_groups.items():
        if len(gsamples) > 1:
            batches = set([])
            for d in gsamples:
                batches |= set(dd.get_batches(d))
            gsamples.sort(key=dd.get_variantcaller_order)
            cur = copy.deepcopy(gsamples[0])
            cur.update({"batch_id": sorted(list(batches))[0] if batches else "_".join(batch_samples),
                        "batch_samples": batch_samples,
                        "variants": {"variantcallers": [dd.get_variantcaller(d) for d in gsamples],
                                     "calls": [d.get("vrn_file") for d in gsamples]}})
            out.append(cur)

    def by_original_order(d):
        return min([sample_order.index(s) for s in d["batch_samples"] if s in sample_order])
    return sorted(out, key=by_original_order)

def combine_calls(*args):
    """Combine multiple callsets into a final set of merged calls.
    """
    if len(args) == 3:
        is_cwl = False
        batch_id, samples, data = args
        caller_names, vrn_files = _organize_variants(samples, batch_id)
    else:
        is_cwl = True
        samples = [utils.to_single_data(x) for x in args]
        samples = [cwlutils.unpack_tarballs(x, x) for x in samples]
        data = samples[0]
        batch_id = data["batch_id"]
        caller_names = data["variants"]["variantcallers"]
        vrn_files = data["variants"]["calls"]
    logger.info("Ensemble consensus calls for {0}: {1}".format(
        batch_id, ",".join(caller_names)))
    edata = copy.deepcopy(data)
    base_dir = utils.safe_makedir(os.path.join(edata["dirs"]["work"], "ensemble", batch_id))
    if any([vcfutils.vcf_has_variants(f) for f in vrn_files]):
        # Decompose multiallelic variants and normalize
        passonly = not tz.get_in(["config", "algorithm", "ensemble", "use_filtered"], edata, False)
        vrn_files = [normalize.normalize(f, data, passonly=passonly, rerun_effects=False, remove_oldeffects=True,
                                         nonrefonly=True,
                                         work_dir=utils.safe_makedir(os.path.join(base_dir, c)))
                     for c, f in zip(caller_names, vrn_files)]
        if "classifiers" not in (dd.get_ensemble(edata) or {}):
            callinfo = _run_ensemble_intersection(batch_id, vrn_files, caller_names, base_dir, edata)
        else:
            config_file = _write_config_file(batch_id, caller_names, base_dir, edata)
            callinfo = _run_ensemble(batch_id, vrn_files, config_file, base_dir,
                                     dd.get_ref_file(edata), edata)
            callinfo["vrn_file"] = vcfutils.bgzip_and_index(callinfo["vrn_file"], data["config"])
        # After decomposing multiallelic variants and normalizing, re-evaluate effects
        ann_ma_file, _ = effects.add_to_vcf(callinfo["vrn_file"], data)
        if ann_ma_file:
            callinfo["vrn_file"] = ann_ma_file

        edata["config"]["algorithm"]["variantcaller"] = "ensemble"
        edata["vrn_file"] = callinfo["vrn_file"]
        edata["ensemble_bed"] = callinfo["bed_file"]
        callinfo["validate"] = validate.compare_to_rm(edata)[0][0].get("validate")
    else:
        out_vcf_file = os.path.join(base_dir, "{0}-ensemble.vcf".format(batch_id))
        vcfutils.write_empty_vcf(out_vcf_file, samples=[dd.get_sample_name(d) for d in samples])
        callinfo = {"variantcaller": "ensemble",
                    "vrn_file": vcfutils.bgzip_and_index(out_vcf_file, data["config"]),
                    "bed_file": None}
    if is_cwl:
        callinfo["batch_samples"] = data["batch_samples"]
        callinfo["batch_id"] = batch_id
        return [{"ensemble": callinfo}]
    else:
        return [[batch_id, callinfo]]

def combine_calls_parallel(samples, run_parallel):
    """Combine calls using batched Ensemble approach.
    """
    batch_groups, extras = _group_by_batches(samples, _has_ensemble)
    out = []
    if batch_groups:
        processed = run_parallel("combine_calls", ((b, xs, xs[0]) for b, xs in batch_groups.items()))
        for batch_id, callinfo in processed:
            for data in batch_groups[batch_id]:
                data["variants"].insert(0, callinfo)
                out.append([data])
    return out + extras

def _has_ensemble(data):
    # for tumour-normal calling, a sample may have "ensemble" for the normal
    # sample configured but there won't be any variant files per se
    variants_to_process = (len(data["variants"]) > 1
                           and any([x.get('vrn_file', None) is not None or x.get('vrn_file_batch', None) is not None
                                    for x in data["variants"]]))
    return variants_to_process and dd.get_ensemble(data)

def _group_by_batches(samples, check_fn):
    """Group calls by batches, processing families together during ensemble calling.
    """
    batch_groups = collections.defaultdict(list)
    extras = []
    for data in [x[0] for x in samples]:
        if check_fn(data):
            batch_groups[multi.get_batch_for_key(data)].append(data)
        else:
            extras.append([data])
    return batch_groups, extras

def _organize_variants(samples, batch_id):
    """Retrieve variant calls for all samples, merging batched samples into single VCF.
    """
    caller_names = [x["variantcaller"] for x in samples[0]["variants"]]
    calls = collections.defaultdict(list)
    for data in samples:
        for vrn in data["variants"]:
            calls[vrn["variantcaller"]].append(vrn["vrn_file"])
    data = samples[0]
    vrn_files = []
    for caller in caller_names:
        fnames = calls[caller]
        if len(fnames) == 1:
            vrn_files.append(fnames[0])
        else:
            vrn_files.append(population.get_multisample_vcf(fnames, batch_id, caller, data))
    return caller_names, vrn_files

def _handle_somatic_ensemble(vrn_file, data):
    """For somatic ensemble, discard normal samples and filtered variants from vcfs.

    Only needed for bcbio.variation based ensemble calling.
    """
    if tz.get_in(["metadata", "phenotype"], data, "").lower().startswith("tumor"):
        vrn_file_temp = vrn_file.replace(".vcf", "_tumorOnly_noFilteredCalls.vcf")
        # Select tumor sample and keep only PASS and . calls
        vrn_file = vcfutils.select_sample(in_file=vrn_file, sample=data["name"][1],
                                          out_file=vrn_file_temp,
                                          config=data["config"], filters="PASS,.")
    return vrn_file

def _bcbio_variation_ensemble(vrn_files, out_file, ref_file, config_file, base_dir, data):
    """Run a variant comparison using the bcbio.variation toolkit, given an input configuration.
    """
    vrn_files = [_handle_somatic_ensemble(v, data) for v in vrn_files]
    tmp_dir = utils.safe_makedir(os.path.join(base_dir, "tmp"))
    resources = config_utils.get_resources("bcbio_variation", data["config"])
    jvm_opts = resources.get("jvm_opts", ["-Xms750m", "-Xmx2g"])
    java_args = ["-Djava.io.tmpdir=%s" % tmp_dir]
    cmd = ["bcbio-variation"] + jvm_opts + java_args + \
          ["variant-ensemble", config_file, ref_file, out_file] + vrn_files
    with utils.chdir(base_dir):
        cmd = "%s %s" % (utils.local_path_export(), " ".join(str(x) for x in cmd))
        do.run(cmd, "Ensemble calling: %s" % os.path.basename(base_dir))

def _run_ensemble(batch_id, vrn_files, config_file, base_dir, ref_file, data):
    """Run an ensemble call using merging and SVM-based approach in bcbio.variation
    """
    out_vcf_file = os.path.join(base_dir, "{0}-ensemble.vcf".format(batch_id))
    out_bed_file = os.path.join(base_dir, "{0}-callregions.bed".format(batch_id))
    work_dir = "%s-work" % os.path.splitext(out_vcf_file)[0]
    if not utils.file_exists(out_vcf_file):
        _bcbio_variation_ensemble(vrn_files, out_vcf_file, ref_file, config_file,
                                  base_dir, data)
        if not utils.file_exists(out_vcf_file):
            base_vcf = glob.glob(os.path.join(work_dir, "prep", "*-cfilter.vcf"))[0]
            utils.symlink_plus(base_vcf, out_vcf_file)
    if not utils.file_exists(out_bed_file):
        multi_beds = glob.glob(os.path.join(work_dir, "prep", "*-multicombine.bed"))
        if len(multi_beds) > 0:
            utils.symlink_plus(multi_beds[0], out_bed_file)
    return {"variantcaller": "ensemble",
            "vrn_file": out_vcf_file,
            "bed_file": out_bed_file if os.path.exists(out_bed_file) else None}

def _write_config_file(batch_id, caller_names, base_dir, data):
    """Write YAML configuration to generate an ensemble set of combined calls.
    """
    config_dir = utils.safe_makedir(os.path.join(base_dir, "config"))
    config_file = os.path.join(config_dir, "{0}-ensemble.yaml".format(batch_id))
    algorithm = data["config"]["algorithm"]
    econfig = {"ensemble": algorithm["ensemble"],
               "names": caller_names,
               "prep-inputs": False}
    intervals = validate.get_analysis_intervals(data, None, base_dir)
    if intervals:
        econfig["intervals"] = os.path.abspath(intervals)
    with open(config_file, "w") as out_handle:
        yaml.safe_dump(econfig, out_handle, allow_unicode=False, default_flow_style=False)
    return config_file

def _get_num_pass(data, n):
    """Calculate the number of samples needed to pass ensemble calling.
    """
    numpass = tz.get_in(["config", "algorithm", "ensemble", "numpass"], data)
    if numpass:
        return int(numpass)
    trusted_pct = tz.get_in(["config", "algorithm", "ensemble", "trusted_pct"], data)
    if trusted_pct:
        return int(math.ceil(float(trusted_pct) * n))
    return 2

def _run_ensemble_intersection(batch_id, vrn_files, callers, base_dir, edata):
    """Run intersection n out of x based ensemble method using bcbio.variation.recall.
    """
    out_vcf_file = os.path.join(base_dir, "{0}-ensemble.vcf.gz".format(batch_id))
    if not utils.file_exists(out_vcf_file):
        num_pass = _get_num_pass(edata, len(vrn_files))
        cmd = [
            config_utils.get_program(
                "bcbio-variation-recall", edata["config"]),
            "ensemble",
            "--cores=%s" % edata["config"]["algorithm"].get("num_cores", 1),
            "--numpass", str(num_pass),
            "--names", ",".join(callers)
        ]
        # Remove filtered calls, do not try to rescue, unless configured
        if not tz.get_in(["config", "algorithm", "ensemble", "use_filtered"], edata):
            cmd += ["--nofiltered"]

        with file_transaction(edata, out_vcf_file) as tx_out_file:
            cmd += [tx_out_file, dd.get_ref_file(edata)] + vrn_files
            cmd = "%s && %s" % (utils.get_java_clprep(), " ".join(str(x) for x in cmd))
            do.run(cmd, "Ensemble intersection calling: %s" % (batch_id))
    in_data = utils.deepish_copy(edata)
    in_data["vrn_file"] = out_vcf_file
    return {"variantcaller": "ensemble",
            "vrn_file": out_vcf_file,
            "bed_file": None}
