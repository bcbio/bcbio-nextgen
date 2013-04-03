"""Perform validation of final calls against known reference materials.

Automates the process of checking pipeline results against known valid calls
to identify discordant variants. This provides a baseline for ensuring the
validity of pipeline updates and algorithm changes.
"""
import os

import yaml

from bcbio import utils
from bcbio.variation import ensemble

def _has_validate(data):
    return data.get("vrn_file") and data["config"]["algorithm"].has_key("validate")

def compare_to_rm(data):
    """Compare final variant calls against reference materials of known calls.
    """
    if _has_validate(data):
        vrn_file = os.path.abspath(data["vrn_file"])
        rm_file = os.path.abspath(data["config"]["algorithm"]["validate"])
        sample = data["name"][-1].replace(" ", "_")
        base_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"],
                                                   "validate", sample,
                                                   data["config"]["algorithm"]["variantcaller"]))
        val_config_file = _create_validate_config_file(vrn_file, rm_file, base_dir, data)
        work_dir = os.path.join(base_dir, "work")
        out = {"summary": os.path.join(work_dir, "validate-summary.csv"),
               "concordant": os.path.join(work_dir, "%s-ref-eval-concordance.vcf" % sample),
               "discordant": os.path.join(work_dir, "%s-eval-ref-discordance-annotate.vcf" % sample)}
        if not utils.file_exists(out["concordant"]):
            ensemble.bcbio_variation_comparison(val_config_file, base_dir, data)
        data["validate"] = out
    return data

def _create_validate_config_file(vrn_file, rm_file, base_dir, data):
    config_dir = utils.safe_makedir(os.path.join(base_dir, "config"))
    config_file = os.path.join(config_dir, "validate.yaml")
    with open(config_file, "w") as out_handle:
        out = _create_validate_config(vrn_file, rm_file, base_dir, data)
        yaml.dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return config_file

def _get_analysis_intervals(data):
    for key in ["callable_regions", "variant_regions"]:
        intervals = data["config"]["algorithm"].get(key)
        if intervals:
            return intervals

def _create_validate_config(vrn_file, rm_file, base_dir, data):
    """Create a bcbio.variation configuration input for validation.
    """
    calls = [{"file": rm_file, "name": "ref", "type": "grading-ref"},
             {"file": vrn_file, "name": "eval"}]
    exp = {"sample": data["name"][-1],
           "ref": data["sam_ref"],
           "align": data["work_bam"],
           "approach": "grade",
           "calls": calls}
    intervals = _get_analysis_intervals(data)
    if intervals:
        exp["intervals"] = intervals
    return {"dir": {"base": base_dir, "out": "work", "prep": "work/prep"},
            "experiments": [exp]}
