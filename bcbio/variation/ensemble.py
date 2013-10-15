"""Ensemble methods that create consensus calls from multiple approaches.

This handles merging calls produced by multiple calling methods or
technologies into a single consolidated callset. Uses the bcbio.variation
toolkit: https://github.com/chapmanb/bcbio.variation
"""
import collections
import copy
import glob
import os

import yaml

from bcbio import utils
from bcbio.log import logger
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.variation import population, validate

def combine_calls(batch_id, samples, data):
    """Combine multiple callsets into a final set of merged calls.
    """
    logger.info("Ensemble consensus calls for {0}: {1}".format(
        batch_id, ",".join(x["variantcaller"] for x in samples[0]["variants"])))
    edata = copy.deepcopy(data)
    base_dir = utils.safe_makedir(os.path.join(edata["dirs"]["work"], "ensemble", batch_id))
    caller_names, vrn_files = _organize_variants(samples, batch_id)
    config_file = _write_config_file(batch_id, caller_names, base_dir, edata)
    callinfo = _run_ensemble(batch_id, vrn_files, config_file, base_dir,
                             edata["sam_ref"], edata["config"])
    edata["config"]["algorithm"]["variantcaller"] = "ensemble"
    edata["vrn_file"] = callinfo["vrn_file"]
    edata["ensemble_bed"] = callinfo["bed_file"]
    callinfo["validate"] = validate.compare_to_rm(edata)[0][0].get("validate")
    return [[batch_id, callinfo]]

def combine_calls_parallel(samples, run_parallel):
    """Combine calls using batched Ensemble approach.
    """
    batch_groups, extras = _group_by_batches(samples, _has_ensemble)
    out = []
    if batch_groups:
        processed = run_parallel("combine_calls", ((b, xs, xs[0]) for b, xs in batch_groups.iteritems()))
        for batch_id, callinfo in processed:
            for data in batch_groups[batch_id]:
                data["variants"].insert(0, callinfo)
                out.append([data])
    return out + extras

def _has_ensemble(data):
    return len(data["variants"]) > 1 and data["config"]["algorithm"].has_key("ensemble")

def _group_by_batches(samples, check_fn):
    """Group calls by batches, processing families together during ensemble calling.
    """
    batch_groups = collections.defaultdict(list)
    extras = []
    for data in [x[0] for x in samples]:
        if check_fn(data):
            batch = data.get("metadata", {}).get("batch")
            if batch:
                batch_groups[batch].append(data)
            else:
                assert data["name"][-1] not in batch_groups
                batch_groups[data["name"][-1]] = [data]
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

def _bcbio_variation_ensemble(vrn_files, out_file, ref_file, config_file, base_dir, config):
    """Run a variant comparison using the bcbio.variation toolkit, given an input configuration.
    """
    tmp_dir = utils.safe_makedir(os.path.join(base_dir, "tmp"))
    bv_jar = config_utils.get_jar("bcbio.variation",
                                  config_utils.get_program("bcbio_variation", config, "dir"))
    resources = config_utils.get_resources("bcbio_variation", config)
    jvm_opts = resources.get("jvm_opts", ["-Xms750m", "-Xmx2g"])
    java_args = ["-Djava.io.tmpdir=%s" % tmp_dir]
    cmd = ["java"] + jvm_opts + java_args + ["-jar", bv_jar, "variant-ensemble", config_file,
                                             ref_file, out_file] + vrn_files
    with utils.chdir(base_dir):
        do.run(cmd, "Ensemble calling: %s" % os.path.basename(base_dir))

def _run_ensemble(batch_id, vrn_files, config_file, base_dir, ref_file, config):
    out_vcf_file = os.path.join(base_dir, "{0}-ensemble.vcf".format(batch_id))
    out_bed_file = os.path.join(base_dir, "{0}-callregions.bed".format(batch_id))
    work_dir = "%s-work" % os.path.splitext(out_vcf_file)[0]
    if not utils.file_exists(out_vcf_file):
        _bcbio_variation_ensemble(vrn_files, out_vcf_file, ref_file, config_file,
                                  base_dir, config)
        if not utils.file_exists(out_vcf_file):
            base_vcf = glob.glob(os.path.join(work_dir, "prep", "*-cfilter.vcf"))[0]
            os.symlink(base_vcf, out_vcf_file)
    if not utils.file_exists(out_bed_file):
        multi_beds = glob.glob(os.path.join(work_dir, "prep", "*-multicombine.bed"))
        if len(multi_beds) > 0:
            os.symlink(multi_beds[0], out_bed_file)
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
    intervals = validate.get_analysis_intervals(data)
    if intervals:
        econfig["intervals"] = os.path.abspath(intervals)
    with open(config_file, "w") as out_handle:
        yaml.dump(econfig, out_handle, allow_unicode=False, default_flow_style=False)
    return config_file
