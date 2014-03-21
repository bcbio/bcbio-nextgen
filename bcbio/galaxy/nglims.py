"""Integration with Galaxy nglims.
"""
import copy
import string

from bcbio.galaxy.api import GalaxyApiAccess
from bcbio.illumina import flowcell

def get_runinfo(galaxy_url, galaxy_apikey, run_folder):
    """Retrieve run information for a processed directory from Galaxy nglims API.
    """
    galaxy_api = GalaxyApiAccess(galaxy_url, galaxy_apikey)
    fc_name, fc_date = flowcell.parse_dirname(run_folder)
    return galaxy_api.run_details(fc_name, fc_date)

def flatten_lane_detail(runinfo):
    """Provide flattened lane information with multiplexed barcodes separated.
    """
    out = []
    for ldetail in runinfo["details"]:
        for i, barcode in enumerate(ldetail.get("multiplex", [{}])):
            cur = copy.deepcopy(ldetail)
            cur["name"] = "%s-%s" % (ldetail["name"], i + 1)
            cur["description"] = barcode.get("name", ldetail["description"])
            cur["bc_index"] = barcode.get("sequence", "")
            cur["project_name"] = clean_name(ldetail["project_name"])
            out.append(cur)
    return out

def clean_name(xs):
    final = []
    safec = "_"
    for x in xs:
        if x not in string.ascii_letters + string.digits:
            if len(final) > 0 and final[-1] != safec:
                final.append(safec)
        else:
            final.append(x)
    if final[-1] == safec:
        final = final[:-1]
    return "".join(final)
