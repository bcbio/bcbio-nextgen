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
from bcbio.pipeline import config_utils
from bcbio.solexa.flowcell import get_flowcell_info

def organize(dirs, config, run_info_yaml):
    """Organize run information from a passed YAML file or the Galaxy API.

    Creates the high level structure used for subsequent processing.
    """
    if run_info_yaml and os.path.exists(run_info_yaml):
        logger.info("Found YAML samplesheet, using %s instead of Galaxy API" % run_info_yaml)
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
    lane_name = "%s_%s_%s" % (item["lane"], fc_date, fc_name)
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

def _check_sample_config(items, in_file):
    """Identify common problems in input sample configuration files.
    """
    logger.info("Checking sample YAML configuration: %s" % in_file)
    _check_for_duplicates(items, "lane", _okay_with_multiplex)
    _check_for_duplicates(items, "description", _okay_with_multiplex)
    _check_for_misplaced(items, "algorithm",
                         ["resources", "metadata", "analysis",
                          "description", "genome_build", "lane", "files"])

def _run_info_from_yaml(fc_dir, run_info_yaml, config):
    """Read run information from a passed YAML file.
    """
    with open(run_info_yaml) as in_handle:
        loaded = yaml.load(in_handle)
    fc_name = None
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
    if fc_name is None:
        fc_name, fc_date = _unique_flowcell_info()
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
