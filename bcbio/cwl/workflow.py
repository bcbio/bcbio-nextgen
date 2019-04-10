"""Prepare bcbio workflows from input YAML files.
"""
import copy
import pprint

import six
import toolz as tz

from bcbio.pipeline import alignment

ALWAYS_AVAILABLE = ["description", "resources"]
STRING_DICT = ["config__algorithm__ensemble"]
FLAT_DICT = ["config__algorithm__variantcaller"]

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
                # flatten inputs that are scattered on entering a sub-workflow
                if step.parallel.endswith("-parallel"):
                    if i == 0:
                        inputs = [_flatten_nested_input(v) for v in inputs]
                        wf_scatter = [get_base_id(v["id"]) for v in inputs]
                    else:
                        inputs = [_flatten_nested_input(v) if get_base_id(v["id"]) in wf_scatter else v
                                  for v in inputs]
                wf_inputs, inputs = _merge_wf_inputs(inputs, wf_inputs, wf_outputs,
                                                     step.internal, wf_step.parallel, nested_inputs)
                outputs, file_vs, std_vs = _get_step_outputs(wf_step, wf_step.outputs, file_vs, std_vs)
                wf_steps.append(("step", wf_step.name, wf_step.parallel, inputs,
                                 outputs, wf_step.image, wf_step.programs, wf_step.disk, wf_step.cores,
                                 wf_step.no_files))
                parallel_ids = _find_split_vs(outputs, wf_step.parallel)
                wf_outputs = _merge_wf_outputs(outputs, wf_outputs, wf_step.parallel)
            yield "wf_start", wf_inputs
            for wf_step in wf_steps:
                yield wf_step
            wf_outputs = [v for v in wf_outputs
                          if v["id"] not in set(["%s" % _get_string_vid(x) for x in step.internal])]
            yield "upload", wf_outputs
            wf_outputs, file_vs, std_vs = _get_step_outputs(step, wf_outputs, file_vs, std_vs)
            yield "wf_finish", step.name, step.parallel, wf_inputs, wf_outputs, wf_scatter
            file_vs = _extract_from_subworkflow(file_vs, step)
            std_vs = _extract_from_subworkflow(std_vs, step)
        elif hasattr(step, "expression"):
            inputs, parallel_ids, nested_inputs = _get_step_inputs(step, file_vs, std_vs, parallel_ids)
            outputs, file_vs, std_vs = _get_step_outputs(step, step.outputs, file_vs, std_vs)
            parallel_ids = _find_split_vs(outputs, step.parallel)
            yield ("expressiontool", step.name, inputs, outputs, step.expression, step.parallel)
        else:
            inputs, parallel_ids, nested_inputs = _get_step_inputs(step, file_vs, std_vs, parallel_ids)
            outputs, file_vs, std_vs = _get_step_outputs(step, step.outputs, file_vs, std_vs)
            parallel_ids = _find_split_vs(outputs, step.parallel)
            yield ("step", step.name, step.parallel, inputs, outputs, step.image, step.programs,
                   step.disk, step.cores, step.no_files)
    yield "upload", [_get_upload_output(x, file_vs) for x in final_outputs]

def _merge_wf_inputs(new, out, wf_outputs, to_ignore, parallel, nested_inputs):
    """Merge inputs for a sub-workflow, adding any not present inputs in out.

    Skips inputs that are internally generated or generated and ignored, keeping
    only as inputs those that we do not generate internally.
    """
    internal_generated_ids = []
    for vignore in to_ignore:
        vignore_id = _get_string_vid(vignore)
        # ignore anything we generate internally, but not those we need to pull in
        # from the external process
        if vignore_id not in [v["id"] for v in wf_outputs]:
            internal_generated_ids.append(vignore_id)
    ignore_ids = set(internal_generated_ids + [v["id"] for v in wf_outputs])
    cur_ids = set([v["id"] for v in out])
    remapped_new = []
    for v in new:
        remapped_v = copy.deepcopy(v)
        outv = copy.deepcopy(v)
        outv["id"] = get_base_id(v["id"])
        outv["source"] = v["id"]
        if outv["id"] not in cur_ids and outv["id"] not in ignore_ids:
            if nested_inputs and v["id"] in nested_inputs:
                outv = _flatten_nested_input(outv)
            out.append(outv)
        if remapped_v["id"] in set([v["source"] for v in out]):
            remapped_v["source"] = get_base_id(remapped_v["id"])
        remapped_new.append(remapped_v)
    return out, remapped_new

def _merge_wf_outputs(new, cur, parallel):
    """Merge outputs for a sub-workflow, replacing variables changed in later steps.

    ignore_ids are those used internally in a sub-workflow but not exposed to subsequent steps
    """
    new_ids = set([])
    out = []
    for v in new:
        outv = {}
        outv["source"] = v["id"]
        outv["id"] = "%s" % get_base_id(v["id"])
        outv["type"] = v["type"]
        if "secondaryFiles" in v:
            outv["secondaryFiles"] = v["secondaryFiles"]
        if tz.get_in(["outputBinding", "secondaryFiles"], v):
            outv["secondaryFiles"] = tz.get_in(["outputBinding", "secondaryFiles"], v)
        new_ids.add(outv["id"])
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
        parts = var["id"].split("/")
        if len(parts) > 1:
            if parts[0] in substep_ids:
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

def is_cwl_record(d):
    """Check if an input is a CWL record, from any level of nesting.
    """
    if isinstance(d, dict):
        if d.get("type") == "record":
            return d
        else:
            recs = list(filter(lambda x: x is not None, [is_cwl_record(v) for v in d.values()]))
            return recs[0] if recs else None
    else:
        return None

def _get_step_inputs(step, file_vs, std_vs, parallel_ids, wf=None):
    """Retrieve inputs for a step from existing variables.

    Potentially nests inputs to deal with merging split variables. If
    we split previously and are merging now, then we only nest those
    coming from the split process.
    """
    inputs = []
    skip_inputs = set([])
    for orig_input in [_get_variable(x, file_vs) for x in _handle_special_inputs(step.inputs, file_vs)]:
        inputs.append(orig_input)
    # Only add description and other information for non-record inputs, otherwise batched with records
    if not any(is_cwl_record(x) for x in inputs):
        inputs += [v for v in std_vs if get_base_id(v["id"]) not in skip_inputs]
    nested_inputs = []
    if step.parallel in ["single-merge", "batch-merge"]:
        if parallel_ids:
            inputs = [_nest_variable(x) if x["id"] in parallel_ids else x for x in inputs]
            nested_inputs = parallel_ids[:]
            parallel_ids = []
    elif step.parallel in ["multi-combined"]:
        assert len(parallel_ids) == 0
        nested_inputs = [x["id"] for x in inputs]
        inputs = [_nest_variable(x) for x in inputs]
    elif step.parallel in ["multi-batch"]:
        assert len(parallel_ids) == 0
        nested_inputs = [x["id"] for x in inputs]
        # If we're batching,with mixed records/inputs avoid double nesting records
        inputs = [_nest_variable(x, check_records=(len(inputs) > 1)) for x in inputs]
    # avoid inputs/outputs with the same name
    outputs = [_get_string_vid(x["id"]) for x in step.outputs]
    final_inputs = []
    for input in inputs:
        input["wf_duplicate"] = get_base_id(input["id"]) in outputs
        final_inputs.append(input)
    return inputs, parallel_ids, nested_inputs

def _get_step_outputs(step, outputs, file_vs, std_vs):
    if len(outputs) == 1 and is_cwl_record(outputs[0]):
        # For workflows, no need to create -- inherit from previous step
        if hasattr(step, "workflow"):
            out_rec = copy.deepcopy(outputs[0])
            out_rec["id"] = "%s/%s" % (step.name, out_rec["id"])
            file_output = [out_rec]
        else:
            file_output = [_create_record(outputs[0]["id"], outputs[0].get("fields", []),
                                          step.name, step.inputs, step.unlist, file_vs, std_vs, step.parallel)]
        std_output = []
    else:
        file_output, std_output = _split_variables([_create_variable(x, step, file_vs) for x in outputs])
    if step.parallel in ["single-split", "batch-split", "multi-combined", "multi-batch"]:
        file_output = [_nest_variable(x) for x in file_output]
        std_output = [_nest_variable(x) for x in std_output]
    # For combined output, ensure at the same level (not nested) as main samples
    if step.parallel in ["multi-combined"]:
        file_vs = _merge_variables([_clean_output(_flatten_nested_input(v) if not is_cwl_record(v) else v)
                                    for v in file_output], file_vs)
    else:
        file_vs = _merge_variables([_clean_output(v) for v in file_output], file_vs)
    std_vs = _merge_variables([_clean_output(v) for v in std_output], std_vs)
    return file_output + std_output, file_vs, std_vs

def _flatten_nested_input(v):
    """Flatten a parallel scatter input -- we only get one of them to tools.
    """
    v = copy.deepcopy(v)
    if isinstance(v["type"], dict) and v["type"]["type"] == "array":
        v["type"] = v["type"]["items"]
    else:
        assert isinstance(v["type"], (list, tuple)), v
        new_type = None
        want_null = False
        for x in v["type"]:
            if isinstance(x, dict) and x["type"] == "array":
                new_type = x["items"]
            elif isinstance(x, six.string_types) and x == "null":
                want_null = True
            else:
                new_type = x
        if want_null:
            if not isinstance(new_type, (list, tuple)):
                new_type = [new_type] if new_type is not None else []
            for toadd in ["null", "string"]:
                if toadd not in new_type:
                    new_type.append(toadd)
        assert new_type, v
        v["type"] = new_type
    return v

def _nest_variable(v, check_records=False):
    """Nest a variable when moving from scattered back to consolidated.

    check_records -- avoid re-nesting a record input if it comes from a previous
    step and is already nested, don't need to re-array.
    """
    if (check_records and is_cwl_record(v) and len(v["id"].split("/")) > 1 and
         v.get("type", {}).get("type") == "array"):
        return v
    else:
        v = copy.deepcopy(v)
        v["type"] = {"type": "array", "items": v["type"]}
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
    if isinstance(vid, six.string_types):
        return vid
    assert isinstance(vid, (list, tuple)), vid
    return "__".join(vid)

def _get_variable(vid, variables):
    """Retrieve an input variable from our existing pool of options.
    """
    if isinstance(vid, six.string_types):
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
    from bcbio import structural
    optional = [["config", "algorithm", "coverage"],
                ["config", "algorithm", "variant_regions"],
                ["config", "algorithm", "sv_regions"],
                ["config", "algorithm", "validate"],
                ["config", "algorithm", "validate_regions"]]
    all_vs = set([get_base_id(v["id"]) for v in variables])
    out = []
    for input in inputs:
        if input == ["reference", "aligner", "indexes"]:
            for v in variables:
                vid = get_base_id(v["id"]).split("__")
                if vid[0] == "reference" and vid[1] in alignment.TOOLS:
                    out.append(vid)
        elif input == ["reference", "snpeff", "genome_build"]:
            found_indexes = False
            for v in variables:
                vid = get_base_id(v["id"]).split("__")
                if vid[0] == "reference" and vid[1] == "snpeff":
                    out.append(vid)
                    found_indexes = True
            assert found_indexes, "Found no snpEff indexes in %s" % [v["id"] for v in variables]
        elif input == ["config", "algorithm", "background", "cnv_reference"]:
            for v in variables:
                vid = get_base_id(v["id"]).split("__")
                if (vid[:4] == ["config", "algorithm", "background", "cnv_reference"] and
                      structural.supports_cnv_reference(vid[4])):
                    out.append(vid)
        elif input in optional:
            if _get_string_vid(input) in all_vs:
                out.append(input)
        else:
            out.append(input)
    return out

def _get_upload_output(vid, variables):
    if isinstance(vid, dict) and "id" in vid:
        parent_v = _get_variable(vid["id"], variables)
        v = copy.deepcopy(vid)
        v["id"] = _get_string_vid(vid["id"])
        v["outputSource"] = parent_v["id"]
    else:
        v = _nest_variable(_get_variable(vid, variables))
        v["outputSource"] = v["id"]
        v["id"] = get_base_id(v["id"])
    v.pop("secondaryFiles", None)
    v["type"].pop("secondaryFiles", None)
    return v

def _create_record(name, field_defs, step_name, inputs, unlist, file_vs, std_vs, parallel):
    """Create an output record by rearranging inputs.

    Batching processes create records that reformat the inputs for
    parallelization.
    """
    if field_defs:
        fields = []
        inherit = []
        inherit_all = False
        inherit_exclude = []
        for fdef in field_defs:
            if not fdef.get("type"):
                if fdef["id"] == "inherit":
                    inherit_all = True
                    inherit_exclude = fdef.get("exclude", [])
                else:
                    inherit.append(fdef["id"])
            else:
                cur = {"name": _get_string_vid(fdef["id"]),
                       "type": fdef["type"]}
                fields.append(_add_secondary_to_rec_field(fdef, cur))
        if inherit_all:
            fields.extend(_infer_record_outputs(inputs, unlist, file_vs, std_vs, parallel, exclude=inherit_exclude))
        elif inherit:
            fields.extend(_infer_record_outputs(inputs, unlist, file_vs, std_vs, parallel, inherit))
    else:
        fields = _infer_record_outputs(inputs, unlist, file_vs, std_vs, parallel)
    out = {"id": "%s/%s" % (step_name, name),
           "type": {"name": name,
                    "type": "record",
                    "fields": fields}}
    if parallel in ["batch-single", "multi-batch"]:
        out = _nest_variable(out)
    return out

def _add_secondary_to_rec_field(orig, cur):
    # CWL does not currently support secondaryFiles in fields
    if orig.get("secondaryFiles"):
        cur["secondaryFiles"] = orig.get("secondaryFiles")
    return cur

def _infer_record_outputs(inputs, unlist, file_vs, std_vs, parallel, to_include=None,
                          exclude=None):
    """Infer the outputs of a record from the original inputs
    """
    fields = []
    unlist = set([_get_string_vid(x) for x in unlist])
    input_vids = set([_get_string_vid(v) for v in _handle_special_inputs(inputs, file_vs)])
    to_include = set([_get_string_vid(x) for x in to_include]) if to_include else None
    to_exclude = tuple(set([_get_string_vid(x) for x in exclude])) if exclude else None
    added = set([])
    for raw_v in std_vs + [v for v in file_vs if get_base_id(v["id"]) in input_vids]:
        # unpack record inside this record and un-nested inputs to avoid double nested
        cur_record = is_cwl_record(raw_v)
        if cur_record:
            # unlist = unlist | set([field["name"] for field in cur_record["fields"]])
            nested_vs = [{"id": field["name"], "type": field["type"]} for field in cur_record["fields"]]
        else:
            nested_vs = [raw_v]
        for orig_v in nested_vs:
            if (get_base_id(orig_v["id"]) not in added
                 and (not to_include or get_base_id(orig_v["id"]) in to_include)):
                if to_exclude is None or not get_base_id(orig_v["id"]).startswith(to_exclude):
                    cur_v = {}
                    cur_v["name"] = get_base_id(orig_v["id"])
                    cur_v["type"] = orig_v["type"]
                    if cur_v["name"] in unlist:
                        cur_v = _flatten_nested_input(cur_v)
                    fields.append(_add_secondary_to_rec_field(orig_v, cur_v))
                    added.add(get_base_id(orig_v["id"]))
    return fields

def _create_variable(orig_v, step, variables):
    """Create a new output variable, potentially over-writing existing or creating new.
    """
    # get current variable, and convert to be the output of our process step
    try:
        v = _get_variable(orig_v["id"], variables)
    except ValueError:
        v = copy.deepcopy(orig_v)
        if not isinstance(v["id"], six.string_types):
            v["id"] = _get_string_vid(v["id"])
    for key, val in orig_v.items():
        if key not in ["id", "type"]:
            v[key] = val
    if orig_v.get("type") != "null":
        v["type"] = orig_v["type"]
    v["id"] = "%s/%s" % (step.name, get_base_id(v["id"]))
    return v

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
    """Retrieve the base id for a variant, ignoring workflow/step prefixes.
    """
    return vid.split("/")[-1]

def get_step_prefix(vid):
    parts = vid.split("/")
    if len(parts) > 1:
        return parts[0]

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
            if "items" in cur_type:
                cur_type = cur_type["items"]
            else:
                cur_type = cur_type["type"]
        if (cur_type in ["File", "null", "record"] or
              (isinstance(cur_type, (list, tuple)) and
               ("File" in cur_type or {'items': 'File', 'type': 'array'} in cur_type))):
            file_vs.append(v)
        elif v["id"] in ALWAYS_AVAILABLE:
            std_vs.append(v)
        else:
            file_vs.append(v)
    return file_vs, std_vs
