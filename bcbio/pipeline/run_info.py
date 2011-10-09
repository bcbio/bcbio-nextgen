"""Retrieve run information describing files to process in a pipeline.

This handles two methods of getting processing information: from a Galaxy
next gen LIMS system or an on-file YAML configuration.
"""
import os
import time
import copy
import string
import datetime
import collections

import yaml

from bcbio.pipeline import log
from bcbio.galaxy.api import GalaxyApiAccess
from bcbio.solexa.flowcell import get_flowcell_info

def get_run_info(fc_dir, config, run_info_yaml):
    """Retrieve run information from a passed YAML file or the Galaxy API.
    """
    if run_info_yaml and os.path.exists(run_info_yaml):
        log.info("Found YAML samplesheet, using %s instead of Galaxy API" % run_info_yaml)
        fc_name, fc_date, run_info = _run_info_from_yaml(fc_dir, run_info_yaml)
    else:
        log.info("Fetching run details from Galaxy instance")
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

def _run_info_from_yaml(fc_dir, run_info_yaml):
    """Read run information from a passed YAML file.
    """
    with open(run_info_yaml) as in_handle:
        loaded = yaml.load(in_handle)
    fc_name = None
    try:
        fc_name, fc_date = get_flowcell_info(fc_dir)
    except ValueError:
        pass
    if isinstance(loaded, dict):
        if loaded.has_key("fc_name") and loaded.has_key("fc_date"):
            fc_name = loaded["fc_name"].replace(" ", "_")
            fc_date = str(loaded["fc_date"]).replace(" ", "_")
        loaded = loaded["details"]
    if fc_name is None:
        fc_name, fc_date = _unique_flowcell_info()
    run_details = []
    for i, item in enumerate(loaded):
        if not item.has_key("lane"):
            item["lane"] = _generate_lane(item["files"], i)
        if not item.has_key("description"):
            item["description"] = str(item["lane"])
        run_details.append(item)
    run_info = dict(details=run_details, run_id="")
    return fc_name, fc_date, run_info

def _clean_extra_whitespace(s):
    while s.endswith(("_", "-", " ", ".")):
        s = s[:-1]
    return s

def _generate_lane(fnames, index):
    """Generate a lane identifier from filenames.
    """
    to_remove = ["s_", "sequence"]
    work_names = []
    if isinstance(fnames, str):
        fnames = [fnames]
    for fname in fnames:
        n = os.path.splitext(os.path.basename(fname))[0]
        for r in to_remove:
            n = n.replace(r, "")
        work_names.append(n)
    if len(work_names) == 1:
        return _clean_extra_whitespace(work_names[0])
    else:
        prefix = _clean_extra_whitespace(os.path.commonprefix(work_names))
        return prefix if prefix else str(index+1)

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
