"""Create Common Workflow Language (CWL) runnable files and tools from a world object.
"""
import copy
import json
import operator
import os

import toolz as tz
import yaml

from bcbio import utils

def from_world(world, run_info_file):
    base = utils.splitext_plus(os.path.basename(run_info_file))[0]
    out_dir = utils.safe_makedir("%s-workflow" % (base))
    out_file = os.path.join(out_dir, "%s-main.cwl" % (base))
    samples = [xs[0] for xs in world]  # unpack world data objects
    analyses = list(set([x["analysis"] for x in samples]))
    assert len(analyses) == 1, "Currently support writing CWL for a single analysis type"
    if analyses[0].startswith("variant"):
        prep_variant_cwl(samples, out_dir, out_file)
    else:
        raise NotImplementedError("Unsupported CWL analysis type: %s" % analyses[0])

def _standard_bcbio_cwl(samples, inputs):
    """Retrieve CWL inputs shared amongst different workflows.

    ToDo: calculate system requirements (cores/memory) from input configuration.
    """
    return {"class": "Workflow",
            "hints": [{"class": "DockerRequirement",
                       "dockerImport": "bcbio/bcbio",
                       "dockerImageId": "bcbio/bcbio"}],
            "requirements": [{"class": "EnvVarRequirement",
                              "envDef": [{"envName": "MPLCONFIGDIR", "envValue": "."}]},
                             {"class": "ScatterFeatureRequirement"}],
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
        assert inp["type"]["type"] == "array", inp
        inp_tool = copy.deepcopy(inp)
        inp_tool["inputBinding"] = {"prefix": "%s=" % inp["id"].replace("#", ""), "separate": False,
                                     "itemSeparator": ";;", "position": i}
        inp_tool["type"] = inp["type"]["items"]
        out["inputs"].append(inp_tool)
    # XXX Need to generalize outputs, just a hack for now to test align_prep
    for outp in outputs:
        outp_tool = copy.deepcopy(outp)
        outp_tool["type"] = {"type": "array", "items": "File"}
        outp_tool["default"] = None
        outp_tool["outputBinding"] = {"glob": "align_prep/*.gz",
                                      "secondaryFiles": [".gbi"]}
        out["outputs"].append(outp_tool)
    with open(out_file, "w") as out_handle:
        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return os.path.join("steps", os.path.basename(out_file))

def _step_template(name, step_dir, inputs, outputs, source=""):
    """Templating function for writing a step to avoid repeating namespaces.
    """
    outputs = [x for x in inputs if x["id"].endswith(tuple(outputs))]
    step_file = _write_tool(step_dir, name, inputs, outputs)
    inputs = [{"id": "#%s.%s" % (name, inp["id"].replace("#", "")),
               "source": "#" + ("%s." % source if source else "") + inp["id"].replace("#", "")}
              for inp in inputs]
    return {"run": {"import": step_file},
            "id": "#%s" % name,
            "scatterMethod": "dotproduct",
            "scatter": [x["id"] for x in inputs],
            "inputs": inputs,
            "outputs": [{"id": "#%s.%s" % (name, output["id"].replace("#", ""))} for output in outputs]}

def prep_variant_cwl(samples, out_dir, out_file):
    """Output a CWL decription for running a variant calling workflow.
    """
    step_dir = utils.safe_makedir(os.path.join(out_dir, "steps"))
    sample_json, variables = _flatten_samples(samples, out_file)
    out = _standard_bcbio_cwl(samples, variables)
    s1 = _step_template("prep_align_inputs", step_dir, variables, ["files"])
    out["steps"] = [s1]
    with open(out_file, "w") as out_handle:
        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out_file, sample_json

def _flatten_samples(samples, base_file):
    """Create a flattened JSON representation of data from the bcbio world map.
    """
    out_file = "%s-samples.json" % utils.splitext_plus(base_file)[0]
    flat_data = []
    for data in samples:
        cur_flat = {}
        for key_path in [["description"], ["rgnames"], ["config", "algorithm"], ["metadata"], ["genome_build"],
                         ["files"], ["reference"], ["genome_resources"], ["vrn_file"]]:
            cur_key = "__".join(key_path)
            for flat_key, flat_val in _to_cwldata(cur_key, tz.get_in(key_path, data)):
                cur_flat[flat_key] = flat_val
        flat_data.append(cur_flat)
    out = {}
    for key in sorted(list(set(reduce(operator.add, [d.keys() for d in flat_data])))):
        out[key] = []
        for cur_flat in flat_data:
            out[key].append(cur_flat.get(key))
    with open(out_file, "w") as out_handle:
        json.dump(out, out_handle, sort_keys=True, indent=4, separators=(',', ': '))
        return out_file, _samplejson_to_inputs(out)

def _add_avro_type(inp, val):
    inp["type"] = _get_avro_type(val)
    return inp

def _get_avro_type(val):
    """Infer avro type for the current input.
    """
    if isinstance(val, dict):
        assert val.get("class") == "File"
        return "File"
    elif isinstance(val, (tuple, list)):
        types = []
        for ctype in [_get_avro_type(v) for v in val]:
            if isinstance(ctype, dict):
                nested_types = [x["items"] for x in types if isinstance(x, dict)]
                if ctype["items"] not in nested_types:
                    types.append(ctype)
            elif ctype not in types:
                types.append(ctype)
        return {"type": "array", "items": (types[0] if len(types) == 1 else types)}
    elif val is None:
        return "null"
    else:
        return "string"

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