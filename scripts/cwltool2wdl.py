#!/bin/env python
"""Exploratory code to convert bcbio generated CWL into WDL.

Uses cwltool parser to parse input CWL then, calls out to cwl2wdl
for generation of WDL using this fork of cwl2wdl:

https://github.com/chapmanb/cwl2wdl

Current status:
- Connected workflow, sub-workflow and tool output
- Records -> Object
- Dotproduct scatter parallelization

Needed for WDL finalization to support bcbio:

- structs to define attributes of an Object
  https://github.com/broadinstitute/cromwell/issues/2283
- Implement standard library function to dump input records in a
  JSON format ('write_struct')
- Implement standard library function to read output of
  a struct/object with multiple levels of nesting (`read_struct`)
- Associate secondary files (like `bai`, `tbi`) with primary file ('bam`, 'vcf.gz`)
  https://github.com/broadinstitute/cromwell/issues/2269
"""
from __future__ import print_function
import collections
import os
import subprocess
import sys

from cwl2wdl import generators
from cwl2wdl import base_classes as cwl2wdl_classes
import cwltool.load_tool
import cwltool.workflow

def main(wf_file, json_file):
    wf_file = os.path.abspath(wf_file)
    out_dir = os.path.dirname(wf_file).replace("-workflow", "-wdl")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    main_wf = cwltool.load_tool.load_tool(wf_file, cwltool.workflow.defaultMakeTool)
    records = {}
    main_wf_dict, records = _wf_to_dict(main_wf, records)
    main_wf_dict["structs"] = records
    main_wf_class = cwl2wdl_classes.Workflow(main_wf_dict)
    for wf_class in [x.task_definition for x in main_wf_class.subworkflows] + [main_wf_class]:
        wdl_file = os.path.join(out_dir, "%s.wdl" % wf_class.name)
        wdl_doc = generators.WdlWorkflowGenerator(wf_class).generate_wdl()
        with open(wdl_file, "w") as out_handle:
            out_handle.write(wdl_doc)
    _validate(wdl_file)

def _validate(wdl_file):
    """Run validation on the generated WDL output using wdltool.
    """
    start_dir = os.getcwd()
    os.chdir(os.path.dirname(wdl_file))
    print("Validating", wdl_file)
    subprocess.check_call(["wdltool", "validate", wdl_file])
    os.chdir(start_dir)

def _wf_to_dict(wf, records):
    """Parse a workflow into cwl2wdl style dictionaries for base and sub-workflows.
    """
    inputs, outputs, records = _get_wf_inout(wf, records)
    out = {"name": _id_to_name(_clean_id(wf.tool["id"])), "inputs": inputs,
           "outputs": outputs, "steps": [], "subworkflows": [],
           "requirements": []}
    for step in wf.steps:
        is_subworkflow = isinstance(step.embedded_tool, cwltool.workflow.Workflow)
        inputs, outputs, remapped, prescatter = _get_step_inout(step)
        inputs, scatter = _organize_step_scatter(step, inputs, remapped)
        if is_subworkflow:
            wf_def, records = _wf_to_dict(step.embedded_tool, records)
            out["subworkflows"].append({"id": "%s.%s" % (wf_def["name"], wf_def["name"]), "definition": wf_def,
                                        "inputs": inputs, "outputs": outputs, "scatter": scatter,
                                        "prescatter": prescatter})
        else:
            task_def, records = _tool_to_dict(step.embedded_tool, records, remapped)
            out["steps"].append({"task_id": task_def["name"], "task_definition": task_def,
                                 "inputs": inputs, "outputs": outputs, "scatter": scatter,
                                 "prescatter": prescatter})
    return out, records

def _clean_id(x):
    """Replace non-allowed characters in WDL input
    """
    return x.replace("wf-", "").replace("-", "_")

def _get_step_inout(step):
    """Retrieve set of inputs and outputs connecting steps.
    """
    inputs = []
    outputs = []
    prescatter = collections.defaultdict(list)
    remapped = {}
    assert step.outputs_record_schema["type"] == "record"
    output_names = set([])
    for outp in step.outputs_record_schema["fields"]:
        outputs.append({"id": outp["name"]})
        output_names.add(outp["name"])
    assert step.inputs_record_schema["type"] == "record"
    for inp in step.inputs_record_schema["fields"]:
        source = inp["source"].split("#")[-1].replace("/", ".")
        # Check if we're unpacking from a record, and unpack from our object
        if "valueFrom" in inp:
            attr_access = "['%s']" % inp["name"]
            if inp["valueFrom"].find(attr_access) > 0:
                source += ".%s" % inp["name"]
                if isinstance(inp["type"], dict) and isinstance(inp["type"].get("items"), dict):
                    if inp["type"]["items"].get("type") == "array" and "inputBinding" in inp["type"]:
                        source, prescatter = _unpack_object_array(inp, source, prescatter)
        # Avoid clashing input and output names, WDL requires unique
        if inp["name"] in output_names:
            new_name = inp["name"] + "_input"
            remapped[inp["name"]] = new_name
            inp["name"] = new_name
        inputs.append({"id": inp["name"], "value": source})
    return inputs, outputs, remapped, dict(prescatter)

def _unpack_object_array(inp, source, prescatter):
    """Unpack Array[Object] with a scatter for referencing in input calls.

    There is no shorthand syntax for referencing all items in an array, so
    we explicitly unpack them with a scatter.
    """
    raise NotImplementedError("Currently not used with record/struct/object improvements")
    base_rec, attr = source.rsplit(".", 1)
    new_name = "%s_%s_unpack" % (inp["name"], base_rec.replace(".", "_"))
    prescatter[base_rec].append((new_name, attr, _to_variable_type(inp["type"]["items"])))
    return new_name, prescatter

def _organize_step_scatter(step, inputs, remapped):
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
            scatter_key = remapped.get(scatter_key) or scatter_key
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

def _get_wf_inout(wf, records):
    assert wf.inputs_record_schema["type"] == "record"
    assert wf.outputs_record_schema["type"] == "record"
    inputs = []
    outputs = []
    for inp in wf.inputs_record_schema["fields"]:
        cur_inp, records = _input_to_dict(inp, records)
        inputs.append(cur_inp)
    for outp in wf.outputs_record_schema["fields"]:
        cur_outp, records = _output_to_dict(outp, records)
        outputs.append(cur_outp)
    return inputs, outputs, records

def _record_to_struct(cur_rec, records):
    """Convert a CWL record into a WDL struct/Object.

    Work in progress to support changes to WDL Objects to
    be defined in structs.
    """
    def to_camel_case(x):
        def uppercase_word(w):
            return w[0].upper() + w[1:]
        return "".join(uppercase_word(w) for w in x.split("_"))
    struct_name = to_camel_case(cur_rec["name"].split("/")[-1])
    if struct_name not in records:
        records[struct_name] = collections.OrderedDict()
        for field in cur_rec["fields"]:
            field_type, records = _to_variable_type(field["type"], records)
            records[struct_name][field["name"].split("/")[-1]] = field_type
    return struct_name, records

def _to_variable_type(x, records):
    """Convert CWL variables to WDL variables, handling nested arrays.
    """
    var_mapping = {"string": "String", "File": "File", "null": "String",
                   "long": "Float", "double": "Float", "int": "Int"}
    if isinstance(x, dict):
        if x["type"] == "record":
            struct_name, records = _record_to_struct(x, records)
            return struct_name, records
        else:
            assert x["type"] == "array", x
            cur_type, records = _to_variable_type(x["items"], records)
            return "Array[%s]" % cur_type, records
    elif isinstance(x, (list, tuple)):
        vars = [v for v in x if v != "null"]
        return var_mapping[vars[0]], records
    else:
        return var_mapping[x], records

def _variable_type_to_read_fn(vartype, records):
    """Convert variant types into corresponding WDL standard library functions.
    """
    fn_map = {"String": "read_string", "Array[String]": "read_lines",
              "Array[Array[String]]": "read_tsv",
              "Object": "read_object", "Array[Object]": "read_objects",
              "Array[Array[Object]]": "read_objects",
              "Int": "read_int", "Float": "read_float"}
    for rec_name in records.keys():
        fn_map["%s" % rec_name] = "read_struct"
        fn_map["Array[%s]" % rec_name] = "read_struct"
        fn_map["Array[Array[%s]]" % rec_name] = "read_struct"
    # Read in Files as Strings
    vartype = vartype.replace("File", "String")
    # Can't read arrays of Ints/Floats
    vartype = vartype.replace("Array[Int]", "Array[String]")
    vartype = vartype.replace("Array[Float]", "Array[String]")
    return fn_map[vartype]

def _arg_to_dict(x, requirements):
    if isinstance(x, basestring):
        return {"prefix": "", "position": None, "value": x}
    elif isinstance(x, dict) and "valueFrom" in x and x["valueFrom"].startswith("sentinel_runtime"):
        for r in requirements:
            if r["requirement_type"] == "cpu":
                cores = r["value"]
            elif r["requirement_type"] == "memory":
                ram = "".join(r["value"].split())
        return {"prefix": "", "position": None, "value": "sentinel_runtime=cores,%s,ram,%s" % (cores, ram)}
    else:
        raise NotImplementedError(x)

def _input_to_dict(i, records, remapped=None):
    """Convert CWL input into dictionary required for a cwl2wdl Input object.
    """
    if not remapped: remapped = {}
    var_type, records = _to_variable_type(i["type"], records)
    if var_type.startswith("Array") and "inputBinding" in i.get("type", {}):
        ib = i["type"]["inputBinding"]
    elif "inputBinding" in i:
        ib = i["inputBinding"]
    else:
        ib = {"prefix": None, "itemSeparator": ";;", "position": None}
    name = _id_to_localname(i["id"]) if "id" in i else i["name"]
    return {"name": remapped.get(name) or name,
            "variable_type": var_type,
            "prefix": ib["prefix"], "separator": ib["itemSeparator"],
            "position": ib["position"], "is_required": True,
            "default": i.get("default", None), "separate": ib.get("separate", True)}, records

def _output_to_dict(o, records):
    if "outputSource" in o:
        name = o["outputSource"].split("#")[-1].replace("/", ".")
    elif "id" in o:
        name = _id_to_localname(o["id"])
    else:
        name = o["name"]
    out_file = "wdl.output.%s.txt" % name
    vartype, records = _to_variable_type(o["type"], records)
    read_fn_name = _variable_type_to_read_fn(vartype, records)
    return {"name": name, "variable_type": vartype,
            "output": "%s('%s')" % (read_fn_name, out_file), "is_required": True}, records

def _id_to_localname(input_id):
    return os.path.basename(input_id).split("#")[1]

def _id_to_name(input_id):
    return os.path.splitext(os.path.basename(input_id))[0]

def _tool_to_dict(tool, records, remapped):
    """Parse a tool definition into a cwl2wdl style dictionary.
    """
    requirements = _requirements_to_dict(tool.requirements + tool.hints)
    inputs = []
    outputs = []
    for inp in tool.tool["inputs"]:
        ready_inp, records = _input_to_dict(inp, records, remapped)
        inputs.append(ready_inp)
    for outp in tool.tool["outputs"]:
        ready_outp, records = _output_to_dict(outp, records)
        outputs.append(ready_outp)
    out = {"name": _id_to_name(tool.tool["id"]),
           "baseCommand": " ".join(tool.tool["baseCommand"]),
           "arguments": [_arg_to_dict(a, requirements) for a in tool.tool["arguments"]],
           "inputs": inputs,
           "outputs": outputs,
           "requirements": requirements,
           "stdin": None, "stdout": None}
    return out, records

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
