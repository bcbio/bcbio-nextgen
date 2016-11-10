#!/bin/env python
"""Exploratory code to convert bcbio generated CWL into the Nextflow DSL.

Uses cwltool parser to parse input CWL then converts objects into cwl2wdl
style dictionaries of workflows and processes (https://github.com/chapmanb/cwl2wdl)
before converting into Nextflow.

Current status:
- tools -> processes

ToDo:

- define script inputs from JSON file
- load outputs from tool JSON output
- connected workflows
- sub-workflows
- Batching of inputs/outputs to grouped Channel (or use set)
- Scatter parallelization
"""
from __future__ import print_function
import collections
import os
import pprint
import sys

import cwltool.load_tool
import cwltool.workflow

def main(wf_file, json_file):
    main_wf = cwltool.load_tool.load_tool(wf_file, cwltool.workflow.defaultMakeTool)
    main_wf_dict = _wf_to_dict(main_wf)
    nf_file = "%s.nf" % os.path.splitext(wf_file)[0]
    with open(nf_file, "w") as out_handle:
        for step in main_wf_dict["steps"]:
            nf_step_to_process(step, out_handle)

# --- Nextflow output

_nf_process_tmpl = """
process {task_id} {{
  {directives}
  input:
    {inputs}
  output:
    {outputs}
\"\"\"
{commandline}
\"\"\"
}}
"""

def nf_step_to_process(step, out_handle):
    """Convert CWL step into a nextflow process.
    """
    pprint.pprint(step)

    directives = []
    for req in step["task_definition"]["requirements"]:
        if req["requirement_type"] == "docker":
            directives.append("container '%s'" % req["value"])
        elif req["requirement_type"] == "cpu":
            directives.append("cpus %s" % req["value"])
        elif req["requirement_type"] == "memory":
            directives.append("memory '%s'" % req["value"])

    task_id = step["task_id"]
    directives = "\n  ".join(directives)
    inputs = "\n".join(nf_inputs(step["inputs"], step["task_definition"]["inputs"]))
    outputs = "\n".join(nf_outputs(step["outputs"], step["task_definition"]["outputs"]))
    commandline = (step["task_definition"]["baseCommand"] + " " +
                   " ".join([nf_input_to_cl(i) for i in step["task_definition"]["inputs"]]))
    out_handle.write(_nf_process_tmpl.format(**locals()))

def nf_inputs(inputs, def_inputs):
    return []

def nf_outputs(outputs, def_outputs):
    return []

def nf_input_to_cl(inp):
    """Convert an input description into command line argument.
    """
    sep = " " if inp.get("separate") else ""
    val = "'%s'" % inp.get("default") if inp.get("default") else "$%s" % inp["name"]
    return "%s%s%s" % (inp["prefix"], sep, val)

# --- Parsing CWL to standard representation

def _wf_to_dict(wf):
    """Parse a workflow into cwl2wdl style dictionary.
    """
    inputs, outputs = _get_wf_inout(wf)
    out = {"name": _id_to_name(wf.tool["id"]).replace("-", "_"), "inputs": inputs,
           "outputs": outputs, "steps": [], "subworkflows": [],
           "requirements": []}
    for step in wf.steps:
        inputs, outputs = _get_step_inout(step)
        inputs, scatter = _organize_step_scatter(step, inputs)
        if isinstance(step.embedded_tool, cwltool.workflow.Workflow):
            wf_def = _wf_to_dict(step.embedded_tool)
            out["subworkflows"].append({"id": wf_def["name"], "definition": wf_def,
                                        "inputs": inputs, "outputs": outputs, "scatter": scatter})
        else:
            task_def = _tool_to_dict(step.embedded_tool)
            out["steps"].append({"task_id": task_def["name"], "task_definition": task_def,
                                 "inputs": inputs, "outputs": outputs, "scatter": scatter})
    return out

def _get_step_inout(step):
    """Retrieve set of inputs and outputs connecting steps.
    """
    inputs = []
    outputs = []
    assert step.inputs_record_schema["type"] == "record"
    for inp in step.inputs_record_schema["fields"]:
        source = inp["source"].split("#")[-1].replace("/", ".")
        # Check if we're unpacking from a record, and unpack from our object
        if "valueFrom" in inp:
            attr_access = "['%s']" % inp["name"]
            if inp["valueFrom"].find(attr_access) > 0:
                source += ".%s" % inp["name"]
        inputs.append({"id": inp["name"], "value": source})
    assert step.outputs_record_schema["type"] == "record"
    for outp in step.outputs_record_schema["fields"]:
        outputs.append({"id": outp["name"]})
    return inputs, outputs

def _organize_step_scatter(step, inputs):
    """Add scattering information from inputs, remapping input variables.
    """
    def extract_scatter_id(inp):
        _, ns_var = inp.split("#")
        _, var = ns_var.split("/")
        return var
    scatter_local = {}
    if "scatter" in step.tool:
        assert step.tool["scatterMethod"] == "dotproduct", \
            "Only support dotproduct scattering in conversion to WDL"
        inp_val = collections.OrderedDict()
        for x in inputs:
            inp_val[x["id"]] = x["value"]
        for scatter_key in [extract_scatter_id(x) for x in step.tool["scatter"]]:
            val = inp_val[scatter_key]
            if len(val.split(".")) in [1, 2]:
                base_key = val
                attr = None
            elif len(val.split(".")) == 3:
                orig_location, record, attr = val.split(".")
                base_key = "%s.%s" % (orig_location, record)
            else:
                raise ValueError("Unexpected scatter input: %s" % val)
            local_ref = base_key.split(".")[-1] + "_local"
            scatter_local[base_key] = local_ref
            if attr:
                local_ref += ".%s" % attr
            inp_val[scatter_key] = local_ref
            inputs = [{"id": iid, "value": ival} for iid, ival in inp_val.items()]
    return inputs, [(v, k) for k, v in scatter_local.items()]

def _get_wf_inout(wf):
    assert wf.inputs_record_schema["type"] == "record"
    assert wf.outputs_record_schema["type"] == "record"
    inputs = [_input_to_dict(inp) for inp in wf.inputs_record_schema["fields"]]
    outputs = [_output_to_dict(outp) for outp in wf.outputs_record_schema["fields"]]
    return inputs, outputs

def _to_variable_type(x):
    """Convert CWL variables to WDL variables, handling nested arrays.
    """
    var_mapping = {"string": "String", "File": "File", "null": "String",
                   "long": "Float", "int": "Int"}
    if isinstance(x, dict):
        if x["type"] == "record":
            return "Object"
        else:
            assert x["type"] == "array", x
            return "Array[%s]" % _to_variable_type(x["items"])
    elif isinstance(x, (list, tuple)):
        vars = [v for v in x if v != "null"]
        assert len(vars) == 1, vars
        return var_mapping[vars[0]]
    else:
        return var_mapping[x]

def _input_to_dict(i):
    """Convert CWL input into dictionary required for a cwl2wdl Input object.
    """
    var_type = _to_variable_type(i["type"])
    if var_type.startswith("Array") and "inputBinding" in i.get("type", {}):
        ib = i["type"]["inputBinding"]
    elif "inputBinding" in i:
        ib = i["inputBinding"]
    else:
        ib = {"prefix": None, "itemSeparator": None, "position": None}
    return {"name": _id_to_localname(i["id"]) if "id" in i else i["name"],
            "variable_type": var_type,
            "prefix": ib["prefix"], "separator": ib["itemSeparator"],
            "position": ib["position"], "is_required": True,
            "default": i.get("default", None), "separate": ib.get("separate", True)}

def _output_to_dict(o):
    if "outputSource" in o:
        name = o["outputSource"].split("#")[-1].replace("/", ".")
    elif "id" in o:
        name = _id_to_localname(o["id"])
    else:
        name = o["name"]
    return {"name": name,
            "variable_type": _to_variable_type(o["type"]),
            "output": "read_cwl_json('cwl.output.json', '%s')" % name, "is_required": True}

def _id_to_localname(input_id):
    return os.path.basename(input_id).split("#")[1]

def _id_to_name(input_id):
    return os.path.splitext(os.path.basename(input_id))[0]

def _tool_to_dict(tool):
    """Parse a tool definition into a cwl2wdl style dictionary.
    """
    out = {"name": _id_to_name(tool.tool["id"]),
           "baseCommand": " ".join(tool.tool["baseCommand"]),
           "arguments": [],
           "inputs": [_input_to_dict(i) for i in tool.tool["inputs"]],
           "outputs": [_output_to_dict(o) for o in tool.tool["outputs"]],
           "requirements": _requirements_to_dict(tool.requirements + tool.hints),
           "stdin": None, "stdout": None}
    return out

def _requirements_to_dict(rs):
    """Convert supported requirements into dictionary for output.
    """
    out = []
    added = set([])
    for r in rs:
        if r["class"] == "DockerRequirement" and "docker" not in added:
            added.add("docker")
            out.append({"requirement_type": "docker", "value": r["dockerImageId"]})
        elif r["class"] == "ResourceRequirement":
            if "coresMin" in r and "cpu" not in added:
                added.add("cpu")
                out.append({"requirement_type": "cpu", "value": r["coresMin"]})
            if "ramMin" in r and "memory" not in added:
                added.add("memory")
                out.append({"requirement_type": "memory", "value": "%s MB" % r["ramMin"]})
            if "tmpdirMin" in r and "disks" not in added:
                added.add("disks")
                out.append({"requirement_type": "disks", "value": "local-disk %s HDD" % r["tmpdirMin"]})
    return out

if __name__ == "__main__":
    main(*sys.argv[1:])
