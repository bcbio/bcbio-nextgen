"""Ensemble methods that create consensus calls from multiple approaches.

This handles merging calls produced by multiple calling methods or
technologies into a single consolidated callset. Uses the bcbio.variation
toolkit: https://github.com/chapmanb/bcbio.variation
"""
import os
import glob
import copy
import subprocess

import yaml

from bcbio import utils
from bcbio.log import logger
from bcbio.pipeline import config_utils

def combine_calls(data):
    """Combine multiple callsets into a final set of merged calls.
    """
    if len(data["variants"]) > 1 and data["config"]["algorithm"].has_key("ensemble"):
        logger.info("Ensemble consensus calls for {0}: {1}".format(
            ",".join(x["variantcaller"] for x in data["variants"]), data["work_bam"]))
        sample = data["name"][-1].replace(" ", "_")
        base_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "ensemble"))
        config_file = _write_config_file(data, sample, base_dir, "ensemble")
        callinfo = _run_bcbio_variation(config_file, base_dir, sample, data)
        data = copy.deepcopy(data)
        data["variants"].insert(0, callinfo)
        _write_config_file(data, sample, base_dir, "compare")
    return [[data]]

def _run_bcbio_variation(config_file, base_dir, sample, data):
    tmp_dir = utils.safe_makedir(os.path.join(base_dir, "tmp"))
    out_vcf_file = os.path.join(base_dir, "{0}-ensemble.vcf".format(sample))
    out_bed_file = os.path.join(base_dir, "{0}-callregions.bed".format(sample))
    if not utils.file_exists(out_vcf_file):
        bv_jar = config_utils.get_jar("bcbio.variation",
                                      config_utils.get_program("bcbio_variation",
                                                               data["config"], "dir"))
        java_args = ["-Djava.io.tmpdir=%s" % tmp_dir]
        subprocess.check_call(["java"] + java_args + ["-jar", bv_jar, "variant-compare", config_file])
        base_vcf = glob.glob(os.path.join(base_dir, sample, "work", "prep",
                                          "*-cfilter.vcf"))[0]
        base_bed = glob.glob(os.path.join(base_dir, sample, "work", "prep",
                                          "*-multicombine.bed"))[0]
        os.symlink(base_vcf, out_vcf_file)
        os.symlink(base_bed, out_bed_file)

    return {"variantcaller": "ensemble",
            "vrn_file": out_vcf_file,
            "bed_file": out_bed_file}

def _write_config_file(data, sample, base_dir, config_name):
    """Write YAML configuration to generate an ensemble set of combined calls.
    """
    sample_dir = os.path.join(base_dir, sample)
    config_dir = utils.safe_makedir(os.path.join(sample_dir, "config"))
    config_file = os.path.join(config_dir, "{0}.yaml".format(config_name))
    prep_fns = {"ensemble": _prep_config_ensemble, "compare": _prep_config_compare}

    econfig = prep_fns[config_name](sample, data["variants"],
                                    data["work_bam"], data["sam_ref"], sample_dir,
                                    data["config"]["algorithm"].get("variant_regions", None),
                                    data["config"]["algorithm"])
    with open(config_file, "w") as out_handle:
        yaml.dump(econfig, out_handle, allow_unicode=False, default_flow_style=False)
    return config_file

def _prep_config_compare(sample, variants, align_bam, ref_file, base_dir,
                         intervals, algorithm):
    """Write YAML bcbio.variation configuration input for results comparison.

    Preps a config file making it easy to compare finalized combined calls
    to individual inputs.
    """
    return _prep_config_shared(sample, variants, align_bam, ref_file, base_dir,
                               intervals, algorithm, "compare", False)

def _prep_config_ensemble(sample, variants, align_bam, ref_file, base_dir,
                          intervals, algorithm):
    """Prepare a YAML configuration file describing the sample inputs.
    """
    return _prep_config_shared(sample, variants, align_bam, ref_file, base_dir,
                               intervals, algorithm, "work", True)

def _prep_config_shared(sample, variants, align_bam, ref_file, base_dir,
                          intervals, algorithm, work_dir, do_combo):
    combo_name = "combo"
    exp = {"sample": sample, "ref": ref_file, "align": align_bam, "calls": []}
    if do_combo:
        cparams = algorithm["ensemble"].get("classifier-params", {})
        exp["finalize"] = \
          [{"method": "multiple",
            "target": combo_name},
            {"method": "recal-filter",
             "target": [combo_name, variants[0]["variantcaller"]],
             "params": {"support": combo_name,
                        "classifiers": algorithm["ensemble"]["classifiers"],
                        "classifier-type": cparams.get("type", "svm"),
                        "normalize": cparams.get("normalize", "default"),
                        "log-attrs": cparams.get("log-attrs", []),
                        "xspecific": True,
                        "trusted":
                        {"total": algorithm["ensemble"].get("trusted-pct", 0.65)}}}]
    if intervals:
        exp["intervals"] = os.path.abspath(intervals)
    for i, v in enumerate(variants):
        cur = {"name": v["variantcaller"], "file": v["vrn_file"],
               "remove-refcalls": True}
        if algorithm.get("ploidy", 2) == 1:
            cur["make-haploid"] = True
        # add a recall variant for the first sample which will combine all calls
        if i == 0 and do_combo:
            recall = copy.deepcopy(cur)
            recall["name"] = combo_name
            recall["recall"] = True
            recall["annotate"] = True
            if algorithm["ensemble"].get("format-filters"):
                recall["format-filters"] = algorithm["ensemble"]["format-filters"]
            exp["calls"].append(recall)
        exp["calls"].append(cur)
    return {"dir": {"base": base_dir, "out": work_dir, "prep": os.path.join(work_dir, "prep")},
            "experiments": [exp]}
