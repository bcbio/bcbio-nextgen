"""Retrieve run information describing files to process in a pipeline.

This handles two methods of getting processing information: from a Galaxy
next gen LIMS system or an on-file YAML configuration.
"""
import copy
import datetime
import itertools
import os
import string
import time

import yaml

from bcbio.log import logger
from bcbio.galaxy.api import GalaxyApiAccess
from bcbio.pipeline import alignment, config_utils
from bcbio.solexa.flowcell import get_flowcell_info
from bcbio.variation import genotype

def organize(dirs, config, run_info_yaml):
    """Organize run information from a passed YAML file or the Galaxy API.

    Creates the high level structure used for subsequent processing.
    """
    if run_info_yaml and os.path.exists(run_info_yaml):
        logger.info("Using input YAML configuration: %s" % run_info_yaml)
        run_details = _run_info_from_yaml(dirs["flowcell"], run_info_yaml, config)
    else:
        logger.info("Fetching run details from Galaxy instance")
        fc_name, fc_date = get_flowcell_info(dirs["flowcell"])
        galaxy_api = GalaxyApiAccess(config['galaxy_url'], config['galaxy_api_key'])
        run_details = []
        galaxy_info = galaxy_api.run_details(fc_name, fc_date)
        for item in galaxy_info["details"]:
            item["upload"] = {"method": "galaxy", "run_id": galaxy_info["run_id"],
                              "fc_name": fc_name, "fc_date": fc_date}
            run_details.append(item)
    out = []
    for item in run_details:
        item["config"] = config_utils.update_w_custom(config, item)
        item["dirs"] = dirs
        if "name" not in item:
            item["name"] = ["", item["description"]]
        out.append(item)
    return out

def _clean_characters(x):
    """Clean problem characters in sample lane or descriptions.
    """
    for problem in [" "]:
        x = x.replace(problem, "_")
    return x

def prep_rg_names(item, config, fc_name, fc_date):
    """Generate read group names from item inputs.
    """
    if fc_name and fc_date:
        lane_name = "%s_%s_%s" % (item["lane"], fc_date, fc_name)
    else:
        lane_name = item["description"]
    return {"rg": item["lane"],
            "sample": item["description"],
            "lane": lane_name,
            "pl": item.get("algorithm", {}).get("platform",
                    config.get("algorithm", {}).get("platform", "illumina")).lower(),
            "pu": lane_name}

def _check_for_duplicates(xs, attr, check_fn=None):
    """Identify and raise errors on duplicate items.
    """
    dups = []
    for key, vals in itertools.groupby(x[attr] for x in xs):
        if len(list(vals)) > 1:
            dups.append(key)
    if len(dups) > 0:
        psamples = []
        for x in xs:
            if x[attr] in dups:
                psamples.append(x)
        # option to skip problem based on custom input function.
        if check_fn and check_fn(psamples):
            return
        descrs = [x["description"] for x in psamples]
        raise ValueError("Duplicate '%s' found in input sample configuration.\n"
                         "Required to be unique for a project: %s\n"
                         "Problem found in these samples: %s" % (attr, dups, descrs))

def _okay_with_multiplex(xs):
    for x in xs:
        if "multiplex" not in x and "barcode" not in x:
            return False
    return True

def _check_for_misplaced(xs, subkey, other_keys):
    """Ensure configuration keys are not incorrectly nested under other keys.
    """
    problems = []
    for x in xs:
        check_dict = x.get(subkey, {})
        for to_check in other_keys:
            if to_check in check_dict:
                problems.append((x["description"], to_check, subkey))
    if len(problems) > 0:
        raise ValueError("\n".join(["Incorrectly nested keys found in sample YAML. These should be top level:",
                                    " sample         |   key name      |   nested under ",
                                    "----------------+-----------------+----------------"] +
                                   ["% 15s | % 15s | % 15s" % (a, b, c) for (a, b, c) in problems]))

def _check_algorithm_keys(item):
    """Check for unexpected keys in the algorithm section.

    Needs to be manually updated when introducing new keys, but avoids silent bugs
    with typos in key names.
    """
    url = "https://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#algorithm-parameters"
    supported = set(["platform", "aligner", "bam_clean", "bam_sort", "trim_reads", "adapters",
                     "custom_trim", "align_split_size", "quality_bin", "quality_format",
                     "write_summary", "merge_bamprep", "coverage", "coverage_bigwig",
                     "coverage_depth", "coverage_interval", "hybrid_target", "hybrid_bait",
                     "ploidy",
                     "variantcaller", "variant_regions", "mark_duplicates",
                     "recalibrate", "realign",
                     "phasing", "validate", "validate_regions", "validate_genome_build",
                     "clinical_reporting",
                     "nomap_split_size", "nomap_split_targets",
                     "ensemble"])
    problem_keys = [k for k in item["algorithm"].iterkeys() if k not in supported]
    if len(problem_keys) > 0:
        raise ValueError("Unexpected configuration keyword in 'algorithm' section: %s\n"
                         "See configuration documentation for supported options:\n%s\n"
                         % (problem_keys, url))

def _check_aligner(item):
    """Ensure specified aligner is valid choice.
    """
    allowed = set(alignment._tools.keys() + [None, False])
    if item["algorithm"].get("aligner") not in allowed:
        raise ValueError("Unexpected algorithm 'aligner' parameter: %s\n"
                         "Supported options: %s\n" %
                         (item["algorithm"].get("aligner"), sorted(list(allowed))))

def _check_variantcaller(item):
    """Ensure specified variantcaller is a valid choice.
    """
    allowed = set(genotype.get_variantcallers().keys() + [None, False])
    vcs = item["algorithm"].get("variantcaller", "gatk")
    if not isinstance(vcs, (tuple, list)):
        vcs = [vcs]
    problem = [x for x in vcs if x not in allowed]
    if len(problem) > 0:
        raise ValueError("Unexpected algorithm 'variantcaller' parameter: %s\n"
                         "Supported options: %s\n" % (problem, sorted(list(allowed))))

def _check_sample_config(items, in_file):
    """Identify common problems in input sample configuration files.
    """
    logger.info("Checking sample YAML configuration: %s" % in_file)
    _check_for_duplicates(items, "lane", _okay_with_multiplex)
    _check_for_duplicates(items, "description", _okay_with_multiplex)
    _check_for_misplaced(items, "algorithm",
                         ["resources", "metadata", "analysis",
                          "description", "genome_build", "lane", "files"])
    [_check_algorithm_keys(x) for x in items]
    [_check_aligner(x) for x in items]
    [_check_variantcaller(x) for x in items]

def _run_info_from_yaml(fc_dir, run_info_yaml, config):
    """Read run information from a passed YAML file.
    """
    with open(run_info_yaml) as in_handle:
        loaded = yaml.load(in_handle)
    fc_name, fc_date = None, None
    if fc_dir:
        try:
            fc_name, fc_date = get_flowcell_info(fc_dir)
        except ValueError:
            pass
    global_config = {}
    if isinstance(loaded, dict):
        global_config = copy.deepcopy(loaded)
        del global_config["details"]
        if loaded.has_key("fc_name") and loaded.has_key("fc_date"):
            fc_name = loaded["fc_name"].replace(" ", "_")
            fc_date = str(loaded["fc_date"]).replace(" ", "_")
        loaded = loaded["details"]
    run_details = []
    for i, item in enumerate(loaded):
        if not item.has_key("lane"):
            item["lane"] = str(i+1)
        item["lane"] = _clean_characters(str(item["lane"]))
        if not item.has_key("description"):
            item["description"] = str(item["lane"])
        item["description"] = _clean_characters(str(item["description"]))
        item["description_filenames"] = global_config.get("description_filenames", False)
        upload = global_config.get("upload", {})
        if fc_name and fc_date:
            upload["fc_name"] = fc_name
            upload["fc_date"] = fc_date
        upload["run_id"] = ""
        item["upload"] = upload
        item["rgnames"] = prep_rg_names(item, config, fc_name, fc_date)
        run_details.append(item)
    _check_sample_config(run_details, run_info_yaml)
    return run_details

def _unique_flowcell_info():
    """Generate data and unique identifier for non-barcoded flowcell.

    String encoding from:
    http://stackoverflow.com/questions/561486/
    how-to-convert-an-integer-to-the-shortest-url-safe-string-in-python
    """
    alphabet = string.ascii_uppercase + string.digits
    fc_date = datetime.datetime.now().strftime("%y%m%d")
    n = int(time.time())
    s = []
    while True:
        n, r = divmod(n, len(alphabet))
        s.append(alphabet[r])
        if n == 0: break
    return ''.join(reversed(s)), fc_date
