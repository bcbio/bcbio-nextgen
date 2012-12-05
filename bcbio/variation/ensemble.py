"""Ensemble methods that create consensus calls from multiple approaches.

This handles merging calls produced by multiple calling methods or
technologies into a single consolidated callset. Uses the bcbio.variation
toolkit: https://github.com/chapmanb/bcbio.variation
"""
import os
import copy
import yaml

from bcbio import utils
from bcbio.log import logger

def combine_calls(data):
    """Combine multiple callsets into a final set of merged calls.
    """
    if len(data["variants"]) > 1 and data["config"]["algorithm"].has_key("ensemble"):
        logger.info("Ensemble consensus calls for {0}: {1}".format(
            ",".join(x["variantcaller"] for x in data["variants"]), data["work_bam"]))
        base_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "ensemble"))
        config_file = _write_config_file(data, base_dir)
        print config_file
        raise NotImplementedError
    return [[data]]

def _write_config_file(data, base_dir):
    config_dir = utils.safe_makedir(os.path.join(base_dir, "config"))
    config_file = os.path.join(config_dir, "ensemble.yaml")
    econfig = _prep_ensemble_config(data["name"][-1], data["variants"],
                                    data["work_bam"], data["sam_ref"], base_dir,
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
                                     "trusted": {"total": algorithm["ensemble"].get("trusted_pct", 0.65)}}}]}
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
