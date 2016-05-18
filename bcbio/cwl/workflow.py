"""Prepare bcbio workflows from input YAML files.
"""
import copy
import pprint

import toolz as tz

ALWAYS_AVAILABLE = ["#description"]

def generate(variables, steps, final_outputs):
    """Generate all of the components of a CWL workflow from input steps.

    file_vs and std_vs are the list of world variables, split into those that
    reference files (and need declaration at each step) and those that don't
    and can be safely passed to each step. This function keeps track of these
    as they're added and updated by steps, passing this information on to
    creation of the CWL files.
    """
    file_vs, std_vs = _split_variables([_flatten_nested_input(v) for v in variables])
    parallel_ids = []
    for step in steps:
        if hasattr(step, "workflow"):
            wf_inputs = []
            wf_outputs = []
            wf_steps = []
            for i, wf_step in enumerate(step.workflow):
                inputs, parallel_ids, nested_inputs = _get_step_inputs(wf_step, file_vs, std_vs, parallel_ids, step)
                outputs, file_vs, std_vs = _get_step_outputs(wf_step, wf_step.outputs, file_vs, std_vs)
                parallel_ids = _find_split_vs(outputs, wf_step.parallel)
                wf_steps.append(("step", wf_step.name, wf_step.parallel, [_clean_record(x) for x in inputs],
                                 outputs, wf_step.programs, wf_step.disk))
                wf_inputs = _merge_wf_inputs(inputs, wf_inputs, wf_outputs, step.internal, wf_step.parallel,
                                             nested_inputs)
                wf_outputs = _merge_wf_outputs(outputs, wf_outputs, wf_step.parallel)
            yield "wf_start", wf_inputs
            for wf_step in wf_steps:
                yield wf_step
            wf_outputs = [v for v in wf_outputs
                          if v["id"] not in set(["#%s" % _get_string_vid(x) for x in step.internal])]
            yield "upload", wf_outputs
            wf_outputs, file_vs, std_vs = _get_step_outputs(step, wf_outputs, file_vs, std_vs)
            yield "wf_finish", step.name, step.parallel, wf_inputs, wf_outputs
            file_vs = _extract_from_subworkflow(file_vs, step)
            std_vs = _extract_from_subworkflow(std_vs, step)
        else:
            inputs, parallel_ids, nested_inputs = _get_step_inputs(step, file_vs, std_vs, parallel_ids)
            outputs, file_vs, std_vs = _get_step_outputs(step, step.outputs, file_vs, std_vs)
            parallel_ids = _find_split_vs(outputs, step.parallel)
            yield "step", step.name, step.parallel, inputs, outputs, step.programs, step.disk
    yield "upload", [_get_upload_output(x, file_vs) for x in final_outputs]

def _merge_wf_inputs(new, out, wf_outputs, to_ignore, parallel, nested_inputs):
    """Merge inputs for a sub-workflow, adding any not present inputs in out.

    Skips inputs that are internally generated or generated and ignored, keeping
    only as inputs those that we do not generate interally.
    """
    internal_generated_ids = []
    for vignore in to_ignore:
        vignore_id = "#%s" % _get_string_vid(vignore)
        # ignore anything we generate internally, but not those we need to pull in
        # from the external process'
        if vignore_id in [v["id"] for v in wf_outputs]:
            internal_generated_ids.append(vignore_id)
    ignore_ids = set(internal_generated_ids + [v["id"] for v in wf_outputs])
    cur_ids = set([v["id"] for v in out])
    for v in new:
        outv = copy.deepcopy(v)
        outv["id"] = "#%s" % get_base_id(v["id"])
        if outv["id"] not in cur_ids and outv["id"] not in ignore_ids:
            if nested_inputs and v["id"] in nested_inputs:
                outv = _flatten_nested_input(outv)
            out.append(outv)
    return out

def _merge_wf_outputs(new, cur, parallel):
    """Merge outputs for a sub-workflow, replacing variables changed in later steps.

    ignore_ids are those used internally in a sub-workflow but not exposed to subsequent steps
    """
    new_ids = set([])
    out = []
    for v in new:
        outv = {}
        outv["source"] = v["id"]
        outv["id"] = "#%s" % get_base_id(v["id"])
        outv["type"] = v["type"]
        if "secondaryFiles" in v:
            outv["secondaryFiles"] = v["secondaryFiles"]
        if tz.get_in(["outputBinding", "secondaryFiles"], v):
            outv["secondaryFiles"] = tz.get_in(["outputBinding", "secondaryFiles"], v)
        new_ids.add(outv["id"])
        if parallel in ["single-split", "batch-split"]:
            outv = _flatten_nested_input(outv)
        out.append(outv)
    for outv in cur:
        if outv["id"] not in new_ids:
            out.append(outv)
    return out

def _extract_from_subworkflow(vs, step):
    """Remove internal variable names when moving from sub-workflow to main.
    """
    substep_ids = set([x.name for x in step.workflow])
    out = []
    for var in vs:
        internal = False
        parts = var["id"].split(".")
        if len(parts) > 1:
            ns = parts[0].replace("#", "")
            if ns in substep_ids:
                internal = True
        if not internal:
            var.pop("source", None)
            out.append(var)
    return out

def _find_split_vs(out_vs, parallel):
    """Find variables created by splitting samples.
    """
    # split parallel job
    if parallel in ["single-parallel", "batch-parallel"]:
        return [v["id"] for v in out_vs]
    else:
        return []

def _get_step_inputs(step, file_vs, std_vs, parallel_ids, wf=None):
    """Retrieve inputs for a step from existing variables.

    Potentially nests inputs to deal with merging split variables. If
    we split previously and are merging now, then we only nest those
    combing from the split process.
    """
    skip_inputs = set([_get_string_vid(x) for x in step.noinputs])
    if wf:
        skip_inputs = skip_inputs | set([_get_string_vid(x) for x in wf.noinputs])
    inputs = []
    for orig_input in [_get_variable(x, file_vs) for x in _handle_special_inputs(step.inputs, file_vs)]:
        is_record_input = tz.get_in(["type", "type"], orig_input) == "record"
        if is_record_input:
            unpack_inputs = _unpack_record(orig_input, skip_inputs)
            inputs.extend(unpack_inputs)
            skip_inputs = skip_inputs | set([get_base_id(v["id"]) for v in unpack_inputs])
        elif get_base_id(orig_input["id"]) not in skip_inputs:
            inputs.append(orig_input)
    inputs += [v for v in std_vs if get_base_id(v["id"]) not in skip_inputs]
    nested_inputs = []
    if step.parallel in ["single-merge", "batch-merge"]:
        if parallel_ids:
            inputs = [_nest_variable(x) if x["id"] in parallel_ids else x for x in inputs]
            nested_inputs = parallel_ids[:]
            parallel_ids = []
    elif step.parallel in ["multi-combined", "multi-batch"]:
        assert len(parallel_ids) == 0
        nested_inputs = [x["id"] for x in inputs]
        inputs = [_nest_variable(x) for x in inputs]
    return inputs, parallel_ids, nested_inputs

def _unpack_record(rec, noinputs):
    """Unpack a record object, extracting individual elements.
    """
    out = []
    for field in rec["type"]["fields"]:
        if field["name"] not in noinputs:
            out.append({"id": "#%s" % field["name"], "type": field["type"],
                        "source": rec["id"], "valueFrom": "$(self.%s)" % field["name"]})
    return out

def _clean_record(var):
    """Remove record source information from an input variant.
    """
    out = copy.deepcopy(var)
    for attr in ["source", "valueFrom"]:
        out.pop(attr, None)
    return out

def _get_step_outputs(step, outputs, file_vs, std_vs):
    if step.parallel in ["multi-batch"]:
        file_output = [_create_record(outputs, step.inputs, step.noinputs, step.unlist, file_vs, std_vs)]
        std_output = []
    else:
        file_output, std_output = _split_variables([_create_variable(x, step, file_vs) for x in outputs])
        #file_output = [_clean_output_extras(x) for x in file_output]
    std_vs = _merge_variables([_clean_output(v) for v in std_output], std_vs)
    file_vs = _merge_variables([_clean_output(v) for v in file_output], file_vs)
    if step.parallel in ["single-split", "batch-split", "multi-combined", "multi-batch"]:
        file_output = [_nest_variable(x) for x in file_output]
        std_output = [_nest_variable(x) for x in std_output]
    return file_output + std_output, file_vs, std_vs

def _flatten_nested_input(v):
    """Flatten a parallel scatter input -- we only get one of them to tools.
    """
    v = copy.deepcopy(v)
    assert v["type"]["type"] == "array"
    v["type"] = v["type"]["items"]
    return v

def _nest_variable(v):
    """Nest a variable when moving from scattered back to consolidated.
    """
    v = copy.deepcopy(v)
    v["type"] = {"type": "array", "items": v["type"]}
    if "secondaryFiles" in v:
        v["type"]["secondaryFiles"] = v.pop("secondaryFiles")
    return v

def _clean_output(v):
    """Remove output specific variables to allow variables to be inputs to next steps.
    """
    out = copy.deepcopy(v)
    outb = out.pop("outputBinding", {})
    if "secondaryFiles" in outb:
        out["secondaryFiles"] = outb["secondaryFiles"]
    return out

def _get_string_vid(vid):
    assert isinstance(vid, (list, tuple)), vid
    return "__".join(vid)

def _get_variable(vid, variables):
    """Retrieve an input variable from our existing pool of options.
    """
    if isinstance(vid, basestring):
        vid = get_base_id(vid)
    else:
        vid = _get_string_vid(vid)
    for v in variables:
        if vid == get_base_id(v["id"]):
            return copy.deepcopy(v)
    raise ValueError("Did not find variable %s in \n%s" % (vid, pprint.pformat(variables)))

def _handle_special_inputs(inputs, variables):
    """Adjust input variables based on special cases.

    This case handles inputs where we are optional or can have flexible choices.

    XXX Need to better expose this at a top level definition.
    """
    optional = [["config", "algorithm", "coverage"],
                ["config", "algorithm", "variant_regions"],
                ["config", "algorithm", "validate"],
                ["config", "algorithm", "validate_regions"]]
    all_vs = set([get_base_id(v["id"]) for v in variables])
    out = []
    for input in inputs:
        if input == ["reference", "aligner", "indexes"]:
            found_indexes = False
            for v in variables:
                vid = get_base_id(v["id"]).split("__")
                if vid[0] == "reference" and vid[-1] == "indexes":
                    out.append(vid)
                    found_indexes = True
            assert found_indexes, "Found no aligner indexes in %s" % [v["id"] for v in variables]
        elif input in optional:
            if _get_string_vid(input) in all_vs:
                out.append(input)
        else:
            out.append(input)
    return out

def _clean_output_extras(v):
    """Remove extra variables we don't want in output.

    We want secondaryFiles nested in outputBinding.
    """
    out = copy.deepcopy(v)
    for not_in_output in ["secondaryFiles"]:
        if not_in_output in out and not_in_output in out.get("outputBinding", {}):
            out.pop(not_in_output, None)
    return out

def _get_upload_output(vid, variables):
    v = _nest_variable(_get_variable(vid, variables))
    v["source"] = v["id"]
    v["id"] = "#%s" % get_base_id(v["id"])
    v["linkMerge"] = "merge_flattened"
    v.pop("secondaryFiles", None)
    return _clean_output_extras(v)

def _create_record(name, inputs, noinputs, unlist, file_vs, std_vs):
    """Create an input record created from rearranging inputs.

    Batching processes create records that reformat the inputs for
    parallelization.
    """
    skip_inputs = set([_get_string_vid(x) for x in noinputs])
    unlist = set([_get_string_vid(x) for x in unlist])
    fields = []
    input_vids = set([_get_string_vid(v) for v in inputs])
    for orig_v in std_vs + [v for v in file_vs if get_base_id(v["id"]) in input_vids]:
        if get_base_id(orig_v["id"]) not in skip_inputs:
            cur_v = {}
            cur_v["name"] = get_base_id(orig_v["id"])
            cur_v["type"] = orig_v["type"]
            if cur_v["name"] in unlist:
                cur_v = _flatten_nested_input(cur_v)
            fields.append(_nest_variable(cur_v))

    return {"id": "#%s" % name,
            "type": {"name": name,
                     "type": "record",
                     "fields": fields}}

def _create_variable(orig_v, step, variables):
    """Create a new output variable, potentially over-writing existing or creating new.
    """
    # get current variable, and convert to be the output of our process step
    try:
        v = _get_variable(orig_v["id"], variables)
    except ValueError:
        v = copy.deepcopy(orig_v)
        if not isinstance(v["id"], basestring):
            v["id"] = "#" + _get_string_vid(v["id"])
    for key, val in orig_v.items():
        if key not in ["id", "type"]:
            v[key] = val
    if v.get("type") == "null" and orig_v.get("type") != "null":
        v["type"] = orig_v["type"]
    return _convert_to_step_id(v, step)

def _merge_variables(new, cur):
    """Add any new variables to the world representation in cur.

    Replaces any variables adjusted by previous steps.
    """
    new_added = set([])
    out = []
    for cur_var in cur:
        updated = False
        for new_var in new:
            if get_base_id(new_var["id"]) == get_base_id(cur_var["id"]):
                out.append(new_var)
                new_added.add(new_var["id"])
                updated = True
                break
        if not updated:
            out.append(cur_var)
    for new_var in new:
        if new_var["id"] not in new_added:
            out.append(new_var)
    return out

def get_base_id(vid):
    """Retrieve the base id for a variant, ignoring prefixes and steps.
    """
    assert vid[0] == "#", vid
    parts = vid[1:].split(".")
    return parts[-1]

def _convert_to_step_id(v, step):
    """Convert a variable to include the identifier of a specific step.
    """
    v["id"] = "#%s.%s" % (step.name, get_base_id(v["id"]))
    return v

def _split_variables(variables):
    """Split variables into always passed (std) and specified (file).

    We always pass some variables to each step but need to
    explicitly define file and algorithm variables so they can
    be linked in as needed.
    """
    file_vs = []
    std_vs = []
    for v in variables:
        cur_type = v["type"]
        while isinstance(cur_type, dict):
            cur_type = cur_type["items"]
        if (cur_type == "File" or cur_type == "null" or
              (isinstance(cur_type, (list, tuple)) and
               ("File" in cur_type or {'items': 'File', 'type': 'array'} in cur_type))):
            file_vs.append(v)
        elif v["id"] in ALWAYS_AVAILABLE:
            std_vs.append(v)
        else:
            file_vs.append(v)
    return file_vs, std_vs
