"""Retrieve run information describing files to process in a pipeline.

This handles two methods of getting processing information: from a Galaxy
next gen LIMS system or an on-file YAML configuration.
"""
import os
import time
import string
import datetime

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
    return fc_name, fc_date, run_info

def _run_info_from_yaml(fc_dir, run_info_yaml):
    """Read run information from a passed YAML file.
    """
    try:
        fc_name, fc_date = get_flowcell_info(fc_dir)
    except ValueError:
        fc_name, fc_date = _unique_flowcell_info()
    with open(run_info_yaml) as in_handle:
        loaded = yaml.load(in_handle)
    run_details = []
    for i, item in enumerate(loaded):
        if not item.has_key("lane"):
            item["lane"] = _generate_lane(item["files"], i)
        run_details.append(item)
    lanes = [x["lane"] for x in run_details]
    assert len(lanes) == len(set(lanes)), "Non unique lanes: %s" % lanes
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
