"""Create Common Workflow Language (CWL) runnable files and tools from a world object.
"""
import json
import os

import toolz as tz
import yaml

from bcbio import utils

def from_world(world, run_info_file):
    base = utils.splitext_plus(os.path.basename(run_info_file))[0]
    out_dir = utils.safe_makedir("%s-workflow" % (base))
    out_file = os.path.join(out_dir, "%s-main.cwl" % (base))
    samples = [xs[0] for xs in world]  # unpack world data objects
    analyses = set([x["analysis"] for x in samples])
    assert len(analyses) == 1, "Currently support writing CWL for a single analysis type"
    if analyses.pop().startswith("variant"):
        prep_variant_cwl(samples, out_dir, out_file)
    else:
        raise NotImplementedError("Unsupported CWL analysis type: %s" % analyses.pop())

def _standard_bcbio_cwl(samples, inputs):
    """Retrieve CWL inputs shared amongst different workflows.

    ToDo: calculate system requirements (cores/memory) from input configuration.
    """
    return {"class": "Workflow",
            "hints": [{"class": "DockerRequirement",
                       "dockerImport": "https://s3.amazonaws.com/bcbio_nextgen/bcbio-nextgen-docker-image.gz",
                       "dockerImageId": "chapmanb/bcbio-nextgen-devel"}],
            "requirements": [{"class": "EnvVarRequirement",
                              "envDef": [{"envName": "MPLCONFIGDIR", "envValue": "."}]}],
            "inputs": inputs,
            "outputs": [],
            "steps": []}

def _write_tool(step_dir, name, inputs, outputs):
    out_file = os.path.join(step_dir, "%s.cwl" % name)
    out = {"class": "CommandLineTool",
           "baseCommand": ["bcbio_nextgen.py", "runfn", name, "cwl"],
           "inputs": [],
           "outputs": []}
    for i, inp in enumerate(inputs):
        out["inputs"].append(tz.assoc(inp, "inputBinding",
                                      {"prefix": "%s=" % inp["id"].replace("#", ""), "separate": False,
                                       "itemSeparator": ";;", "position": i}))
    with open(out_file, "w") as out_handle:
        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return os.path.join("steps", os.path.basename(out_file))

def _step_template(name, step_dir, inputs, outputs, source=""):
    """Templating function for writing a step to avoid repeating namespaces.
    """
    step_file = _write_tool(step_dir, name, inputs, outputs)
    outputs = [x for x in inputs if x["id"].endswith(tuple(outputs))]
    return {"run": {"import": step_file},
            "id": "#%s" % name,
            "inputs": [{"id": "#%s.%s" % (name, inp["id"].replace("#", "")),
                        "source": "#" + ("%s." % source if source else "") + inp["id"].replace("#", "")}
                       for inp in inputs],
            "outputs": [{"id": "#%s.%s" % (name, output["id"].replace("#", ""))} for output in outputs]}

def prep_variant_cwl(samples, out_dir, out_file):
    """Output a CWL decription for running a variant calling workflow.
    """
    step_dir = utils.safe_makedir(os.path.join(out_dir, "steps"))
    sample_json, variables = _flatten_samples(samples, out_file)
    out = _standard_bcbio_cwl(samples, variables)
    s1 = _step_template("prep_align_inputs", step_dir, variables, ["/files"])
    out["steps"] = [s1]
    with open(out_file, "w") as out_handle:
        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out_file, sample_json

def _flatten_samples(samples, base_file):
    """Create a flattened JSON representation of data from the bcbio world map.
    """
    out_file = "%s-samples.json" % utils.splitext_plus(base_file)[0]
    out = {}
    for data in samples:
        ns = data["description"]
        for key_path in [["rgnames"], ["config", "algorithm"], ["metadata"], ["genome_build"], ["files"],
                         ["reference"], ["genome_resources"], ["vrn_file"]]:
            cur_key = "%s__%s" % (ns, "__".join(key_path))
            for flat_key, flat_val in _to_cwldata(cur_key, tz.get_in(key_path, data)):
                out[flat_key] = flat_val
    with open(out_file, "w") as out_handle:
        json.dump(out, out_handle, sort_keys=True, indent=4, separators=(',', ': '))
        return out_file, _samplejson_to_inputs(out)

def _add_avro_type(inp, val):
    """Infer avro type for the current input.
    """
    if isinstance(val, dict):
        assert val.get("class") == "File"
        inp["type"] = "File"
    elif isinstance(val, (tuple, list)):
        if isinstance(val[0], dict):
            assert val[0].get("class") == "File"
            cur_type = "File"
        else:
            cur_type = "string"
        inp["type"] = {"type": "array", "items": cur_type}
    elif val is None:
        inp["type"] = "null"
    else:
        inp["type"] = "string"
    return inp

def _samplejson_to_inputs(svals):
    """Convert sample output into inputs for CWL configuration files, with types.
    """
    out = []
    for key, val in svals.items():
        out.append(_add_avro_type({"id": "#%s" % key}, val))
    return out

def _to_cwldata(key, val):
    """Convert nested dictionary into CWL data, flatening and marking up files.

    Moves file objects to the top level, enabling insertion in CWL inputs/outputs.
    """
    out = []
    if isinstance(val, dict):
        remain_val = {}
        for nkey, nval in val.items():
            cur_nkey = "%s__%s" % (key, nkey)
            cwl_nval = _item_to_cwldata(nval)
            if isinstance(cwl_nval, dict):
                out.extend(_to_cwldata(cur_nkey, nval))
            elif cwl_nval == nval:
                remain_val[nkey] = nval
            else:
                out.append((cur_nkey, cwl_nval))
        if remain_val:
            out.append((key, json.dumps(remain_val, separators=(',', ':'))))
    else:
        out.append((key, _item_to_cwldata(val)))
    return out

def _item_to_cwldata(x):
    """"Markup an item with CWL specific metadata.
    """
    if isinstance(x, (list, tuple)):
        return [_item_to_cwldata(subx) for subx in x]
    elif x and isinstance(x, basestring) and (os.path.isfile(x) or os.path.isdir(x)) and os.path.exists(x):
        return {"class": "File", "path": x}
    else:
        return x