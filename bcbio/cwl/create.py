"""Create Common Workflow Language (CWL) runnable files and tools from a world object.
"""
import copy
import json
import operator
import os

import toolz as tz
import yaml

from bcbio import utils
from bcbio.cwl import workflow

def from_world(world, run_info_file):
    base = utils.splitext_plus(os.path.basename(run_info_file))[0]
    out_dir = utils.safe_makedir("%s-workflow" % (base))
    out_file = os.path.join(out_dir, "main-%s.cwl" % (base))
    samples = [xs[0] for xs in world]  # unpack world data objects
    analyses = list(set([x["analysis"] for x in samples]))
    assert len(analyses) == 1, "Currently support writing CWL for a single analysis type"
    if analyses[0].startswith("variant"):
        prep_variant_cwl(samples, out_dir, out_file)
    else:
        raise NotImplementedError("Unsupported CWL analysis type: %s" % analyses[0])

def _cwl_workflow_template(inputs):
    """Retrieve CWL inputs shared amongst different workflows.
    """
    return {"class": "Workflow",
            "hints": [{"class": "DockerRequirement",
                       "dockerPull": "bcbio/bcbio",
                       "dockerImageId": "bcbio/bcbio"}],
            "requirements": [{"class": "EnvVarRequirement",
                              "envDef": [{"envName": "MPLCONFIGDIR", "envValue": "."}]},
                             {"class": "ScatterFeatureRequirement"},
                             {"class": "SubworkflowFeatureRequirement"},
                             {"class": "InlineJavascriptRequirement"}],
            "inputs": inputs,
            "outputs": [],
            "steps": []}

def _write_tool(step_dir, name, inputs, outputs, parallel):
    out_file = os.path.join(step_dir, "%s.cwl" % name)
    out = {"class": "CommandLineTool",
           "baseCommand": ["bcbio_nextgen.py", "runfn", name, "cwl"],
           "inputs": [],
           "outputs": []}
    pinputs = [{"id": "#sentinel-parallel", "type": "string",
                "default": parallel}]
    inputs = pinputs + inputs
    for i, inp in enumerate(inputs):
        base_id = workflow.get_base_id(inp["id"])
        inp_tool = copy.deepcopy(inp)
        inp_tool["id"] = "#%s" % base_id
        inp_binding = {"prefix": "%s=" % base_id, "separate": False,
                       "itemSeparator": ";;", "position": i}
        inp_tool = _place_input_binding(inp_tool, inp_binding, parallel)
        inp_tool = _place_secondary_files(inp_tool, inp_binding)
        out["inputs"].append(inp_tool)
    for outp in outputs:
        outp_tool = copy.deepcopy(outp)
        outp_tool["id"] = "#%s" % workflow.get_base_id(outp["id"])
        out["outputs"].append(outp_tool)
    with open(out_file, "w") as out_handle:
        def str_presenter(dumper, data):
            if len(data.splitlines()) > 1:  # check for multiline string
                return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')
            return dumper.represent_scalar('tag:yaml.org,2002:str', data)
        yaml.add_representer(str, str_presenter)
        yaml.dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return os.path.join("steps", os.path.basename(out_file))

def _place_input_binding(inp_tool, inp_binding, parallel):
    """Check nesting of variables to determine where to place the input binding.

    We want to allow having multiple files together (like fasta_indices), combined
    with the itemSeparator, but also support having multiple samples where we pass
    things independently.
    """
    if parallel == "multi-combined" and tz.get_in(["type", "type"], inp_tool) == "array":
        inp_tool["type"]["inputBinding"] = inp_binding
    else:
        inp_tool["inputBinding"] = inp_binding
    return inp_tool

def _place_secondary_files(inp_tool, inp_binding):
    """Put secondaryFiles at the level of the File item to ensure indexes get passed.

    This involves using a second input binding to get the secondaryFiles, that
    we ignore downstream. Ideally we could use `valueFrom: null` but that doesn't
    seem to work right now.
    """
    secondary_files = inp_tool.pop("secondaryFiles", None)
    if secondary_files:
        key = []
        while tz.get_in(key + ["type"], inp_tool) != "File" and tz.get_in(key + ["items"], inp_tool) != "File":
            key.append("type")
        secondary_key = key + ["inputBinding"]
        if tz.get_in(secondary_key, inp_tool):
            inp_tool = tz.update_in(inp_tool, secondary_key + ["secondaryFiles"], lambda x: secondary_files)
        else:
            nested_inp_binding = copy.deepcopy(inp_binding)
            nested_inp_binding["prefix"] = "ignore="
            nested_inp_binding["secondaryFiles"] = secondary_files
            inp_tool = tz.update_in(inp_tool, secondary_key, lambda x: nested_inp_binding)
    return inp_tool

def _step_template(name, run_file, inputs, outputs, parallel):
    """Templating function for writing a step to avoid repeating namespaces.
    """
    scatter_inputs = []
    sinputs = []
    for inp in inputs:
        step_inp = {"id": "#%s.%s" % (name, workflow.get_base_id(inp["id"])), "source": inp["id"]}
        sinputs.append(step_inp)
        # scatter on inputs from previous processes that have been arrayed
        if parallel == "multi-parallel" or len(inp["id"].split(".")) > 1:
            scatter_inputs.append(step_inp["id"])
    out = {"run": {"import": run_file},
           "id": "#%s" % name,
           "inputs": sinputs,
           "outputs": [{"id": output["id"]} for output in outputs]}
    if parallel in ["single-parallel", "multi-parallel"]:
        out.update({"scatterMethod": "dotproduct",
                    "scatter": scatter_inputs})
    return out

def prep_variant_cwl(samples, out_dir, out_file):
    """Output a CWL decription for running a variant calling workflow.
    """
    step_dir = utils.safe_makedir(os.path.join(out_dir, "steps"))
    sample_json, variables = _flatten_samples(samples, out_file)
    out = _cwl_workflow_template(variables)
    parent_wfs = []
    for cur in workflow.variant(variables):
        if cur[0] == "step":
            _, name, parallel, inputs, outputs = cur
            step_file = _write_tool(step_dir, name, inputs, outputs, parallel)
            out["steps"].append(_step_template(name, step_file, inputs, outputs, parallel))
        elif cur[0] == "upload":
            out["outputs"] = cur[1]
        elif cur[0] == "wf_start":
            parent_wfs.append(out)
            out = _cwl_workflow_template(cur[1])
        elif cur[0] == "wf_finish":
            _, name, parallel, inputs, outputs = cur
            wf_out_file = "wf-%s.cwl" % name
            with open(os.path.join(out_dir, wf_out_file), "w") as out_handle:
                yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
            out = parent_wfs.pop(-1)
            out["steps"].append(_step_template(name, wf_out_file, inputs, outputs, parallel))
        else:
            raise ValueError("Unexpected workflow value %s" % str(cur))

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
        for key_path in [["analysis"], ["description"], ["rgnames"], ["config", "algorithm"],
                         ["metadata"], ["genome_build"],
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
        assert val.get("class") == "File" or "File" in val.get("class")
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
        out = {"class": "File", "path": x}
        if os.path.isfile(x):
            if x.endswith(".bam"):
                out["secondaryFiles"] = [{"class": "File", "path": x + ".bai"}]
            elif x.endswith((".vcf.gz", ".bed.gz")):
                out["secondaryFiles"] = [{"class": "File", "path": x + ".tbi"}]
        return out
    else:
        return x
