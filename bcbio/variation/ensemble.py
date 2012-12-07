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
from bcbio.pipeline import config_loader as config_utils

def combine_calls(data):
    """Combine multiple callsets into a final set of merged calls.
    """
    if len(data["variants"]) > 1 and data["config"]["algorithm"].has_key("ensemble"):
        logger.info("Ensemble consensus calls for {0}: {1}".format(
            ",".join(x["variantcaller"] for x in data["variants"]), data["work_bam"]))
        sample = data["name"][-1].replace(" ", "_")
        base_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "ensemble"))
        config_file = _write_config_file(data, sample, base_dir)
        callinfo = _run_bcbio_variation(config_file, base_dir, sample, data)
        print callinfo
        data["variants"] = data["variants"].insert(0, callinfo)
    return [[data]]

def _run_bcbio_variation(config_file, base_dir, sample, data):
    out_vcf_file = os.path.join(base_dir, "{0}-ensembl.vcf".format(sample))
    out_bed_file = os.path.join(base_dir, "{0}-callregions.bed".format(sample))
    if not utils.file_exists(out_vcf_file):
        bv_jar = config_utils.get_jar("bcbio.variation",
                                      config_utils.get_program("bcbio_variation",
                                                               data["config"], "dir"))
        subprocess.check_call(["java", "-jar", bv_jar, "variant-compare", config_file])
        base_vcf = glob.glob(os.path.join(base_dir, sample, "work", "prep", "*-cfilter.vcf"))[0]
        base_bed = glob.glob(os.path.join(base_dir, sample, "work", "prep", "*-multicombine.bed"))[0]
        os.symlink(base_vcf, out_vcf_file)
        os.symlink(base_bed, out_bed_file)

    return {"variantcaller": "ensemble",
            "vrn_file": out_vcf_file,
            "bed_file": out_bed_file}

def _write_config_file(data, sample, base_dir):
    sample_dir = os.path.join(base_dir, sample)
    config_dir = utils.safe_makedir(os.path.join(sample_dir, "config"))
    config_file = os.path.join(config_dir, "ensemble.yaml")
    econfig = _prep_ensemble_config(sample, data["variants"],
                                    data["work_bam"], data["sam_ref"], sample_dir,
                                    data["config"]["algorithm"].get("variant_regions", None),
                                    data["config"]["algorithm"])
    with open(config_file, "w") as out_handle:
        yaml.dump(econfig, out_handle, allow_unicode=False, default_flow_style=False)
    return config_file

def _prep_ensemble_config(sample, variants, align_bam, ref_file, base_dir,
                          intervals, algorithm):
    """Prepare a YAML configuration file describing the sample inputs.
    """
    combo_name = "combo"
    exp = {"sample": sample, "ref": ref_file, "align": align_bam, "calls": [],
           "finalize": [{"method": "multiple",
                         "target": combo_name},
                        {"method": "recal-filter",
                         "target": [combo_name, variants[1]["variantcaller"]],
                         "params": {"support": combo_name,
                                    "classifiers": algorithm["ensemble"]["classifiers"],
                                    "xspecific": True,
                                     "trusted":
                                     {"total": algorithm["ensemble"].get("trusted_pct", 0.65)}}}]}
    if intervals:
        exp["intervals"] = os.path.abspath(intervals)
    for i, v in enumerate(variants):
        cur = {"name": v["variantcaller"], "file": v["vrn_file"],
               "annotate": True}
        if algorithm.get("ploidy", 2) == 1:
            cur["make-haploid"] = True
        # add a recall variant for the first sample which will combine all calls
        if i == 0:
            recall = copy.deepcopy(cur)
            recall["name"] = combo_name
            recall["recall"] = True
            if algorithm["ensemble"].get("format-filters"):
                recall["format-filters"] = algorithm["ensemble"]["format-filters"]
            exp["calls"].append(recall)
        exp["calls"].append(cur)
    return {"dir": {"base": base_dir, "out": "work", "prep": "work/prep"},
            "experiments": [exp]}
