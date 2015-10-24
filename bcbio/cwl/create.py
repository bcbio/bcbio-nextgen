"""Create Common Workflow Language (CWL) runnable files and tools from a world object.
"""
import os

import yaml

from bcbio import utils

def from_world(world, run_info_file):
    base = utils.splitext_plus(os.path.basename(run_info_file))[0]
    out_dir = utils.safe_makedir("%s-workflow" % (base))
    out_file = os.path.join(out_dir, "%s-main.cwl" % (base))
    samples = [xs[0] for xs in world] # unpack world data objects
    analyses = set([x["analysis"] for x in samples])
    assert len(analyses) == 1, "Currently support writing CWL for a single analysis type"
    if analyses.pop().startswith("variant"):
        prep_variant_cwl(samples, out_dir, out_file)
    else:
        raise NotImplementedError("Unsupported CWL analysis type: %s" % analyses.pop())

def _standard_bcbio_cwl(samples):
    """Retrieve CWL inputs shared amongst different workflows.

    ToDo: calculate system requirements (cores/memory) from input configuration.
    """
    return {"class": "Workflow",
            "hints": [{"class": "DockerRequirement",
                       "dockerImport": "https://s3.amazonaws.com/bcbio_nextgen/bcbio-nextgen-docker-image.gz",
                       "dockerImageId": "bcbio/bcbio"}],
            "requirements": [{"class": "EnvVarRequirement",
                              "envDef": [{"envName": "MPLCONFIGDIR", "envValue": "."}]}],
            "inputs": [],
            "outputs": [],
            "steps": []}

def _step_template(name, inputs, outputs):
    """Templating function for writing a step to avoid repeating namespaces.
    """
    step_file = "steps/%s.cwl" % name
    return {"run": {"import": step_file},
            "id": "#%s" % name,
            "inputs": [{"id": "#%s.%s" % (name, input), "source": []} for input, source in inputs],
            "outputs": [{"id": "#%s.%s" % (name, output)} for output in outputs]}

def prep_variant_cwl(samples, out_dir, out_file):
    """Output a CWL decription for running a variant calling workflow.
    """
    step_dir = utils.safe_makedir(os.path.join(out_dir, "steps"))
    out = _standard_bcbio_cwl(samples)
    out["steps"] = [_step_template("prep_align_inputs", [("files", ""), ("algorithm", "")], ["files"])]
    with open(out_file, "w") as out_handle:
        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out_file