"""Ensemble methods that create consensus calls from multiple approaches.

This handles merging calls produced by multiple calling methods or
technologies into a single consolidated callset. Uses the bcbio.variation
toolkit: https://github.com/chapmanb/bcbio.variation
"""
import os

from bcbio import utils
from bcbio.log import logger

def combine_calls(data):
    """Combine multiple callsets into a final set of merged calls.
    """
    if len(data["variants"]) > 1 and data["config"]["algorithm"].has_key("ensemble"):
        logger.info("Ensemble consensus calls for {0}: {1}".format(
            ",".join(x["variantcaller"] for x in data["variants"]), data["work_bam"]))
        base_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "ensemble"))
        econfig = _prep_ensemble_config(data["name"][-1], data["variants"],
                                        data["work_bam"], data["sam_ref"], base_dir,
                                        data["config"]["algorithm"].get("variant_regions", None),
                                        data["config"]["algorithm"]["ensemble"])
        raise NotImplementedError
    return [[data]]

def _prep_ensemble_config(sample, variants, align_bam, ref_file, base_dir,
                          intervals, config):
    """Prepare a YAML configuration file describing the sample inputs.
    """
    exp = {"sample": sample, "ref": ref_file, "align": align_bam, "calls": []}
    if intervals:
        exp["intervals"] = intervals
    return {"dir": {"base": base_dir, "out": "work", "prep": "work/prep"},
            "experiments": [exp]}
    return out
