"""Create Common Workflow Language (CWL) runnable files and tools from a world object.
"""
import os

import yaml

from bcbio import utils

def from_world(world_yaml, run_info_file):
    base = utils.splitext_plus(os.path.basename(run_info_file))[0]
    out_dir = utils.safe_makedir("%s-workflow" % (base))
    out_file = os.path.join(out_dir, "%s-main.cwl" % (base))

    out = {"class": "Workflow",
           "requirements": [],
           "inputs": [],
           "outputs": [],
           "steps": []}
    with open(out_file, "w") as out_handle:
        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out_file
