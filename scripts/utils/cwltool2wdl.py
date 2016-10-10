#!/bin/env python
"""Exploratory code to convert bcbio generated CWL into WDL.

Uses cwltool parser to parse input CWL then, calls out to cwl2wdl
for generation of WDL using this fork of cwl2wdl:

https://github.com/chapmanb/cwl2wdl

Targets for supporting in conversion:

- Workflows
- Sub-workflows
- Tools -> Tasks
- scatter parallelization
- Record -> Object/Map

Current status:
- Basic workflow support at top level
- Initial Tool output

ToDo:
- Map inputs in workflows, associating tasks
- Read WDL output for each variable from JSON files (like read_json, but one output file?)
- Figure out how to convert Records to Object/Map
- Convert scatter parallelization
- Add validation of WDL outputs
"""
import os
import pprint
import sys

from cwl2wdl import generators
from cwl2wdl import base_classes as cwl2wdl_classes
import cwltool.load_tool
import cwltool.workflow
import wdl.parser

def main(wf_file, json_file):
    main_wf = cwltool.load_tool.load_tool(wf_file, cwltool.workflow.defaultMakeTool)
    main_wf_dict = cwl2wdl_classes.Workflow(_wf_to_dict(main_wf))
    wdl_doc = generators.WdlWorkflowGenerator(main_wf_dict).generate_wdl()
    wdl_file = "%s.wdl" % os.path.splitext(wf_file)[0]
    with open(wdl_file, "w") as out_handle:
        out_handle.write(wdl_doc)
    #print(wdl.parser.parse(wdl_doc, os.path.basename(wdl_file)))

def _wf_to_dict(wf):
    """Parse a workflow into cwl2wdl style dictionary.
    """
    out = {"name": _id_to_name(wf.tool["id"]), "inputs": [], "outputs": [],
           "steps": [], "subworkflows": [],
           "requirements": []}
    for step in wf.steps:
        if isinstance(step.embedded_tool, cwltool.workflow.Workflow):
            wf_def = _wf_to_dict(step.embedded_tool)
            out["subworkflows"].append({"id": wf_def["name"], "definition": wf_def,
                                        "inputs": [], "outputs": []})
        else:
            task_def = _tool_to_dict(step.embedded_tool)
            out["steps"].append({"task_id": task_def["name"], "task_definition": task_def,
                                 "inputs": [], "outputs": []})
    return out

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
    else:
        ib = i["inputBinding"]
    return {"name": _id_to_localname(i["id"]) , "variable_type": var_type,
            "prefix": ib["prefix"], "separator": ib["itemSeparator"],
            "position": ib["position"], "is_required": True,
            "default": i.get("default", None), "separate": ib.get("separate", True)}

def _id_to_localname(input_id):
    return os.path.basename(input_id).split("#")[1]

def _id_to_name(input_id):
    return os.path.splitext(os.path.basename(input_id))[0]

def _tool_to_dict(tool):
    """Parse a tool definition into a cwl2wdl style dictionary.
    """
    #print(tool)
    #print(dir(tool))
    #pprint.pprint(tool.tool)
    # print(tool.requirements)
    # print(tool.hints)
    out = {"name": _id_to_name(tool.tool["id"]),
           "baseCommand": " ".join(tool.tool["baseCommand"]),
           "arguments": [],
           "inputs": [_input_to_dict(i) for i in tool.tool["inputs"]],
           "outputs": [{"name": _id_to_localname(o["id"]), "variable_type": _to_variable_type(o["type"]),
                        "output": "", "is_required": True}
                       for o in tool.tool["outputs"]],
           "requirements": [],
           "stdin": None, "stdout": None}
    return out

if __name__ == "__main__":
    main(*sys.argv[1:])
