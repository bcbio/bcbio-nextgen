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
    """
    file_vs, std_vs = _split_variables([_flatten_nested_input(v) for v in variables])
    s = collections.namedtuple("s", "name inputs outputs")
    steps = [s("prep_align_inputs", [["files"]],
               [_cwl_file_world(["files"], ".gbi"),
                _cwl_nonfile_world(["config", "algorithm", "quality_format"])]),
             s("process_alignment",
               [["files"], ["reference", "fasta", "indexes"], ["reference", "fasta", "base"],
                ["reference", "bwa", "indexes"]],
               [_cwl_file_world(["work_bam"], ".bai"),
                _cwl_file_world(["align_bam"], ".bai"),
                _cwl_file_world(["work_bam-plus", "disc"], ".bai"),
                _cwl_file_glob(["work_bam-plus", "sr"], "'align/' + name + '/' + name + '-sort-sr.bam'",
                               ".bai")]),
             s("prep_samples",
               [["config", "algorithm", "variant_regions"]],
               [_cwl_file_world(["config", "algorithm", "variant_regions"]),
                _cwl_file_world(["config", "algorithm", "variant_regions_merged"])]),
             s("postprocess_alignment",
               [["align_bam"],
                ["reference", "fasta", "base"], ["reference", "fasta", "indexes"],
                ["config", "algorithm", "variant_regions_merged"]],
               [_cwl_nonfile_world(["config", "algorithm", "coverage_interval"]),
                _cwl_file_world(["regions", "callable"]),
                _cwl_file_world(["regions", "sample_callable"]),
                _cwl_file_world(["regions", "nblock"]),
                _cwl_file_world(["regions", "highdepth"], allow_missing=True),
                _cwl_file_world(["regions", "offtarget_stats"])]),
             # TODO -- combine sample regions should group into sample batches
             s("combine_sample_regions",
               [["regions", "callable"], ["regions", "nblock"],
                ["reference", "fasta", "base"], ["reference", "fasta", "indexes"]],
               [_cwl_file_world(["config", "algorithm", "callable_regions"]),
                _cwl_file_world(["config", "algorithm", "non_callable_regions"]),
                _cwl_nonfile_world(["config", "algorithm", "callable_count"], "int")])]
    for step in steps:
        inputs = [_get_variable(x, file_vs) for x in step.inputs] + std_vs
        file_output, std_output = _split_variables([_create_variable(x, step, file_vs) for x in step.outputs])
        std_vs = _merge_variables([_clean_output(v) for v in std_output], std_vs)
        file_vs = _merge_variables([_clean_output(v) for v in file_output], file_vs)
        yield step.name, inputs, file_output + std_output

def _cwl_nonfile_world(key, keytype="string"):
    """Retrieve a non-file value from a key in the bcbio world object.
    """
    converter = ["return val;"]
    return _cwl_get_from_world(key, converter, keytype)

def _cwl_file_world(key, extension="", allow_missing=False):
    """Retrieve a file, or array of files, from a key in the bcbio world object.
    """
    secondary_str = (", 'secondaryFiles': [{'class': 'File', 'path': dir + val + '%s'}]" % extension) if extension else ""
    converter = ["if (val === null)",
                 "  return val;"
                 "else if (typeof val === 'string' || val instanceof String)",
                 "  return [{'path': dir + val, 'class': 'File'%s}];" % secondary_str,
                 "else",
                 "  var vals = val;",
                 "  return vals.map(function(val){return {'path': dir + val, 'class': 'File'%s};});" % secondary_str]
    return _cwl_get_from_world(key, converter, ["File", 'null'] if allow_missing else "File", extension)

def _cwl_get_from_world(key, convert_val, valtype, extension=""):
    """Generic function to retrieve specific results from a bcbio world object.

    The generic javascript provides `dir`, the directory containing the output files
    potentially remapped externally for Docker containers and `val` -- the value of
    the `key` attribute from the bcbio world object.
    """
    keygetter = "".join(["['%s']" % k for k in key])
    getter = ["${",
              "   var dir = self[0].path.replace(/\/[^\/]*$/,'') + '/';",
              "   var val = JSON.parse(self[0].contents)%s;" % keygetter] + \
              ["   %s" % v for v in convert_val] + \
              ["}"]
    out = {"id": key,
           "type": valtype,
           "outputBinding": {"glob": "cwl-*-world.json",
                             "loadContents": True,
                             "outputEval": "\n".join(getter)}}
    if extension:
        out["outputBinding"]["secondaryFiles"] = [extension]
    return out

def _cwl_file_glob(key, file_pattern, extension=""):
    """Retrieve an output CWL file, and extensions, using glob on the filesystem.
    """
    file_glob = ["${",
                 "   var name = JSON.parse(inputs.rgnames)['sample'];",
                 "   return %s;" % file_pattern,
                 "}"]
    out = {"id": key,
           "type": "File",
           "outputBinding": {"glob": "\n".join(file_glob)}}
    if extension:
        out["outputBinding"]["secondaryFiles"] = [extension]
    return out

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
    vid = _get_string_vid(vid)
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
        v = copy.deepcopy(orig_v)
        v["id"] = "#" + _get_string_vid(v["id"])
    for key, val in orig_v.items():
        if key not in ["id", "type"]:
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
        if cur_type == "File" or isinstance(cur_type, (list, tuple)) and "File" in cur_type:
            file_vs.append(v)
        else:
            std_vs.append(v)
    return file_vs, std_vs