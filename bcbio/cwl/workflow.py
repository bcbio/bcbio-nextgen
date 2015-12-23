"""Prepare bcbio workflows from input YAML files.

This organizes the metadata and other information about workflows,
providing the necessary information to translate into CWL. The goal is to
eventually replace pipeline/main.py workflows with generalized
versions of this code.
"""
import copy
import collections
import pprint

def variant(variables):
    """Variant calling workflow implementation in CWL.

    The logic here is general, aside from `steps` and could be re-used for
    non-variable pipelines. This is an attempt to get the approach right before
    supporting everything with CWL.

    file_vs and std_vs are the list of world variables, split into those that
    reference files (and need declaration at each step) and those that don't
    and can be safely passed to each step. This function keeps track of these
    as they're added and updated by steps, passing this information on to
    creation of the CWL files.

    Working notes/?s:
     - How do we handle cases where there might or might not be an output, like
       prep_align_inputs with potential bgzipping?
    """
    file_vs, std_vs = _split_variables([_flatten_nested_input(v) for v in variables])
    s = collections.namedtuple("s", "name inputs outputs")
    steps = [s("prep_align_inputs", [["files"]],
               [{"id": ["files"],
                 "outputBinding": {"glob": "align_prep/*.gz",
                                   "secondaryFiles": [".gbi"]
                                   }}]),
             s("process_alignment",
               [["files"], ["reference", "fasta", "indexes"], ["reference", "fasta", "base"],
                ["reference", "bwa", "indexes"]],
               [])]
    for step in steps:
        inputs = [_get_variable(x, file_vs) for x in step.inputs] + std_vs
        file_output, std_output = _split_variables([_create_variable(x, step, file_vs) for x in step.outputs])
        std_vs = _merge_variables([_clean_output(v) for v in std_output], std_vs)
        file_vs = _merge_variables([_clean_output(v) for v in file_output], file_vs)
        yield step.name, inputs, file_output + std_output

def _flatten_nested_input(v):
    """Flatten a parallel scatterplot input -- we only get one of them to tools.
    """
    v = copy.deepcopy(v)
    assert v["type"]["type"] == "array"
    v["type"] = v["type"]["items"]
    return v

def _clean_output(v):
    """Remove output specific variables to allow variables to be inputs to next steps.
    """
    out = copy.deepcopy(v)
    for key in ["outputBinding"]:
        out.pop(key, None)
    return out

def _get_variable(vid, variables):
    """Retrieve an input variable from our existing pool of options.
    """
    vid = "__".join(vid)
    for v in variables:
        if vid == get_base_id(v["id"]):
            return copy.deepcopy(v)
    raise ValueError("Did not find variable %s in \n%s" % (vid, pprint.pformat(variables)))

def _create_variable(orig_v, step, variables):
    """Create a new output variable, potentially over-writing existing or creating new.
    """
    # get current variable, and convert to be the output of our process step
    try:
        v = _get_variable(orig_v["id"], variables)
    except ValueError:
        print orig_v, step
        raise NotImplementedError("Need to do creation of non-existing variables")
    for key, val in orig_v.items():
        if key != "id":
            v[key] = val
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
    assert vid[0] == "#"
    parts = vid[1:].split(".")
    return parts[-1]

def _convert_to_step_id(v, step):
    """Convert a variable to include the identifier of a specific step.
    """
    v["id"] = "#%s.%s" % (step.name, get_base_id(v["id"]))
    return v

def _split_variables(variables):
    """Split variables passed into file-based and non-file.

    We always pass non-file variables to each step but need to
    explicitly define file variables so they can be linked in.
    """
    file_vs = []
    std_vs = []
    for v in variables:
        cur_type = v["type"]
        while isinstance(cur_type, dict):
            cur_type = cur_type["items"]
        if cur_type == "File":
            file_vs.append(v)
        else:
            std_vs.append(v)
    return file_vs, std_vs