"""Retrieve run information describing files to process in a pipeline.

This handles two methods of getting processing information: from a Galaxy
next gen LIMS system or an on-file YAML configuration.
"""
import collections
import copy
import datetime
import itertools
import os
import string
import time

import yaml

from bcbio.log import logger
from bcbio.galaxy.api import GalaxyApiAccess
from bcbio.solexa.flowcell import get_flowcell_info

def get_run_info(fc_dir, config, run_info_yaml):
    """Retrieve run information from a passed YAML file or the Galaxy API.
    """
    if run_info_yaml and os.path.exists(run_info_yaml):
        logger.info("Found YAML samplesheet, using %s instead of Galaxy API" % run_info_yaml)
        fc_name, fc_date, run_info = _run_info_from_yaml(fc_dir, run_info_yaml, config)
    else:
        logger.info("Fetching run details from Galaxy instance")
        fc_name, fc_date = get_flowcell_info(fc_dir)
        galaxy_api = GalaxyApiAccess(config['galaxy_url'], config['galaxy_api_key'])
        run_info = galaxy_api.run_details(fc_name, fc_date)
    return fc_name, fc_date, _organize_runs_by_lane(run_info)

def _organize_runs_by_lane(run_info):
    """Organize run information collapsing multiplexed items by lane.

    Lane is the unique identifier in a run and used to combine multiple
    run items on a fastq lane, separable by barcodes.
    """
    items = _normalize_barcodes(run_info["details"])
    items_by_lane = collections.defaultdict(list)
    for x in items:
        items_by_lane[x["lane"]].append(x)
    out = []
    for grouped_items in [items_by_lane[x] for x in sorted(items_by_lane.keys())]:
        bcs = [x["barcode_id"] for x in grouped_items]
        assert len(bcs) == len(set(bcs)), "Duplicate barcodes {0} in lane {1}".format(
            bcs, grouped_items[0]["lane"])
        assert len(bcs) == 1 or None not in bcs, "Barcode and non-barcode in lane {0}".format(
            grouped_items[0]["lane"])
        out.append(grouped_items)
    run_info["details"] = out
    return run_info

def _normalize_barcodes(items):
    """Normalize barcode specification methods into individual items.
    """
    split_items = []
    for item in items:
        if item.has_key("multiplex"):
            for multi in item["multiplex"]:
                base = copy.deepcopy(item)
                base["description"] += ": {0}".format(multi["name"])
                del multi["name"]
                del base["multiplex"]
                base.update(multi)
                split_items.append(base)
        elif item.has_key("barcode"):
            item.update(item["barcode"])
            del item["barcode"]
            split_items.append(item)
        else:
            item["barcode_id"] = None
            split_items.append(item)
    return split_items

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

def _check_sample_config(items):
    """Identify common problems in input sample configuration files.
    """
    _check_for_duplicates(items, "lane", _okay_with_multiplex)
    _check_for_duplicates(items, "description", _okay_with_multiplex)

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
        upload = global_config.get("upload")
        if upload:
            upload["fc_name"] = fc_name
            upload["fc_date"] = fc_date
        item["upload"] = upload
        item["rgnames"] = prep_rg_names(item, config, fc_name, fc_date)
        run_details.append(item)
    _check_sample_config(run_details)
    run_info = dict(details=run_details, run_id="")
    return fc_name, fc_date, run_info

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
