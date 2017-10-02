"""Create Common Workflow Language (CWL) runnable files and tools from a world object.
"""
import copy
import functools
import json
import math
import operator
import os
import tarfile

import toolz as tz
import yaml

from bcbio import utils
from bcbio.cwl import defs, workflow
from bcbio.distributed import objectstore, resources

INTEGRATION_MAP = {"keep:": "arvados", "s3:": "s3", "sbg:": "sbgenomics",
                   "dx:": "dnanexus"}

def from_world(world, run_info_file, integrations=None):
    base = utils.splitext_plus(os.path.basename(run_info_file))[0]
    out_dir = utils.safe_makedir("%s-workflow" % (base))
    out_file = os.path.join(out_dir, "main-%s.cwl" % (base))
    samples = [xs[0] for xs in world]  # unpack world data objects
    analyses = list(set([x["analysis"] for x in samples]))
    assert len(analyses) == 1, "Only support writing CWL for a single analysis type: %s" % analyses
    try:
        workflow_fn = defs.workflows[analyses[0].lower()]
    except KeyError:
        raise NotImplementedError("Unsupported CWL analysis type: %s" % analyses[0])
    prep_cwl(samples, workflow_fn, out_dir, out_file, integrations)

def _cwl_workflow_template(inputs, top_level=False):
    """Retrieve CWL inputs shared amongst different workflows.
    """
    ready_inputs = []
    for inp in inputs:
        cur_inp = copy.deepcopy(inp)
        for attr in ["source", "valueFrom", "wf_duplicate"]:
            cur_inp.pop(attr, None)
        if top_level:
            cur_inp = workflow._flatten_nested_input(cur_inp)
        cur_inp = _clean_record(cur_inp)
        ready_inputs.append(cur_inp)
    return {"class": "Workflow",
            "cwlVersion": "v1.0",
            "hints": [],
            "requirements": [{"class": "EnvVarRequirement",
                              "envDef": [{"envName": "MPLCONFIGDIR", "envValue": "."}]},
                             {"class": "ScatterFeatureRequirement"},
                             {"class": "SubworkflowFeatureRequirement"}],
            "inputs": ready_inputs,
            "outputs": [],
            "steps": []}

def _add_disk_estimates(cwl_res, inputs, file_estimates, disk):
    """Add disk usage estimates to CWL ResourceRequirement.

    Based on inputs (which need to be staged) and disk
    specifications (which estimate outputs).
    """
    if not disk:
        disk = {}
    if file_estimates:
        total_estimate = 0
        for key, multiplier in disk.items():
            if key in file_estimates:
                total_estimate += int(multiplier * file_estimates[key])
        for inp in inputs:
            scale = 2.0 if inp.get("type") == "array" else 1.0
            if workflow.is_cwl_record(inp):
                for f in _get_record_fields(inp):
                    if f["name"] in file_estimates:
                        total_estimate += file_estimates[f["name"]] * scale
            elif inp["id"] in file_estimates:
                total_estimate += file_estimates[inp["id"]] * scale
        if total_estimate:
            # scale total estimate to allow extra room, round to integer
            total_estimate = int(math.ceil(total_estimate * 1.5))
            cwl_res["tmpdirMin"] = total_estimate
            cwl_res["outdirMin"] += total_estimate
    return cwl_res

def _write_tool(step_dir, name, inputs, outputs, parallel, image, programs,
                file_estimates, disk, step_cores, samples):
    out_file = os.path.join(step_dir, "%s.cwl" % name)
    resource_cores, mem_gb_per_core = resources.cpu_and_memory((programs or []) + ["default"], samples)
    cores = min([step_cores, resource_cores]) if step_cores else resource_cores
    mem_mb_total = int(mem_gb_per_core * cores * 1024)
    bcbio_docker_disk = 1 * 1024  # Minimum requirements for bcbio Docker image
    cwl_res = {"class": "ResourceRequirement",
               "coresMin": cores, "ramMin": mem_mb_total, "outdirMin": bcbio_docker_disk}
    cwl_res = _add_disk_estimates(cwl_res, inputs, file_estimates, disk)
    docker_image = "bcbio/bcbio" if image == "bcbio" else "quay.io/bcbio/%s" % image
    docker = {"class": "DockerRequirement", "dockerPull": docker_image, "dockerImageId": docker_image}
    out = {"class": "CommandLineTool",
           "cwlVersion": "v1.0",
           "baseCommand": ["bcbio_nextgen.py", "runfn", name, "cwl"],
           "requirements": [],
           "hints": [docker, cwl_res],
           "arguments": [],
           "inputs": [],
           "outputs": []}
    if programs:
        def resolve_package(p):
            out = {}
            parts = p.split("=")
            if len(parts) == 2:
                out["package"] = parts[0]
                out["version"] = [parts[1]]
            else:
                out["package"] = p
            out["specs"] = ["https://anaconda.org/bioconda/%s" % out["package"]]
            return out
        out["hints"].append({"class": "SoftwareRequirement",
                             "packages": [resolve_package(p) for p in programs]})
    # Use JSON for inputs, rather than command line arguments
    # Correctly handles multiple values and batching across CWL runners
    use_commandline_args = False
    out["requirements"] += [{"class": "InlineJavascriptRequirement"},
                            {"class": "InitialWorkDirRequirement",
                                "listing": [{"entryname": "cwl.inputs.json",
                                            "entry": "$(JSON.stringify(inputs))"}]}]
    out["arguments"] += [{"position": 0, "valueFrom":
                          "sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])"},
                         "sentinel_parallel=%s" % parallel,
                         "sentinel_outputs=%s" % ",".join([_get_sentinel_val(v) for v in outputs]),
                         "sentinel_inputs=%s" % ",".join(["%s:%s" %
                                                          (workflow.get_base_id(v["id"]),
                                                           "record" if workflow.is_cwl_record(v) else "var")
                                                          for v in inputs])]
    for i, inp in enumerate(inputs):
        base_id = workflow.get_base_id(inp["id"])
        inp_tool = copy.deepcopy(inp)
        inp_tool["id"] = base_id
        if inp.get("wf_duplicate"):
            inp_tool["id"] += "_toolinput"
        for attr in ["source", "valueFrom", "wf_duplicate"]:
            inp_tool.pop(attr, None)
        if _is_scatter_parallel(parallel) and _do_scatter_var(inp, parallel):
            inp_tool = workflow._flatten_nested_input(inp_tool)
        if use_commandline_args:
            inp_binding = {"prefix": "%s=" % base_id,
                           "separate": False, "itemSeparator": ";;", "position": i}
            inp_tool = _place_input_binding(inp_tool, inp_binding, parallel)
        else:
            inp_binding = None
        inp_tool = _place_secondary_files(inp_tool, inp_binding)
        inp_tool = _clean_record(inp_tool)
        out["inputs"].append(inp_tool)
    for outp in outputs:
        outp_tool = copy.deepcopy(outp)
        outp_tool = _clean_record(outp_tool)
        outp_tool["id"] = workflow.get_base_id(outp["id"])
        out["outputs"].append(outp_tool)
    with open(out_file, "w") as out_handle:
        def str_presenter(dumper, data):
            if len(data.splitlines()) > 1:  # check for multiline string
                return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')
            return dumper.represent_scalar('tag:yaml.org,2002:str', data)
        yaml.add_representer(str, str_presenter)
        yaml.dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return os.path.join("steps", os.path.basename(out_file))

def _clean_record(rec):
    """Remove secondary files from record fields, which are currently not supported.

    To be removed later when secondaryFiles added to records.
    """
    if workflow.is_cwl_record(rec):
        def _clean_fields(d):
            if isinstance(d, dict):
                if "fields" in d:
                    out = []
                    for f in d["fields"]:
                        f = utils.deepish_copy(f)
                        f.pop("secondaryFiles", None)
                        out.append(f)
                    d["fields"] = out
                    return d
                else:
                    out = {}
                    for k, v in d.items():
                        out[k] = _clean_fields(v)
                    return out
            else:
                return d
        return _clean_fields(rec)
    else:
        return rec

def _get_record_fields(d):
    """Get field names from a potentially nested record.
    """
    if isinstance(d, dict):
        if "fields" in d:
            return d["fields"]
        else:
            for v in d.values():
                fields = _get_record_fields(v)
                if fields:
                    return fields

def _get_sentinel_val(v):
    """Retrieve expected sentinel value for an output, expanding records.
    """
    out = workflow.get_base_id(v["id"])
    if workflow.is_cwl_record(v):
        out += ":%s" % ";".join([x["name"] for x in _get_record_fields(v)])
    return out

def _place_input_binding(inp_tool, inp_binding, parallel):
    """Check nesting of variables to determine where to place the input binding.

    We want to allow having multiple files together (like fasta_indices), combined
    with the itemSeparator, but also support having multiple samples where we pass
    things independently.
    """
    if (parallel in ["multi-combined", "multi-batch", "batch-split", "batch-parallel",
                     "batch-merge", "batch-single"] and
          tz.get_in(["type", "type"], inp_tool) == "array"):
        inp_tool["type"]["inputBinding"] = inp_binding
    else:
        inp_tool["inputBinding"] = inp_binding
    return inp_tool

def _place_secondary_files(inp_tool, inp_binding=None):
    """Put secondaryFiles at the level of the File item to ensure indexes get passed.
    """
    def _is_file(val):
        return (val == "File" or (isinstance(val, (list, tuple)) and "File" in val))
    secondary_files = inp_tool.pop("secondaryFiles", None)
    if secondary_files:
        key = []
        while (not _is_file(tz.get_in(key + ["type"], inp_tool))
               and not _is_file(tz.get_in(key + ["items"], inp_tool))
               and not _is_file(tz.get_in(key + ["items", "items"], inp_tool))):
            key.append("type")
        if tz.get_in(key, inp_tool):
            inp_tool["secondaryFiles"] = secondary_files
        elif inp_binding:
            nested_inp_binding = copy.deepcopy(inp_binding)
            nested_inp_binding["prefix"] = "ignore="
            nested_inp_binding["secondaryFiles"] = secondary_files
            inp_tool = tz.update_in(inp_tool, key, lambda x: nested_inp_binding)
    return inp_tool

def _is_scatter_parallel(parallel):
    return parallel.endswith("-parallel")

def _do_scatter_var(v, parallel):
    """Logic for scattering a variable.
    """
    # For batches, scatter records only at the top level (double nested)
    if parallel.startswith("batch") and workflow.is_cwl_record(v):
        return (tz.get_in(["type", "type"], v) == "array" and
                tz.get_in(["type", "type", "type"], v) == "array")
    # Otherwise, scatter arrays
    else:
        return (tz.get_in(["type", "type"], v) == "array")

def _step_template(name, run_file, inputs, outputs, parallel, scatter=None):
    """Templating function for writing a step to avoid repeating namespaces.
    """
    scatter_inputs = []
    sinputs = []
    for inp in inputs:
        step_inp = {"id": workflow.get_base_id(inp["id"]), "source": inp["id"]}
        if inp.get("wf_duplicate"):
            step_inp["id"] += "_toolinput"
        for attr in ["source", "valueFrom"]:
            if attr in inp:
                step_inp[attr] = inp[attr]
        sinputs.append(step_inp)
        # scatter on inputs from previous processes that have been arrayed
        if (_is_scatter_parallel(parallel) and (_do_scatter_var(inp, parallel)
                                                or (scatter and inp["id"] in scatter))):
            scatter_inputs.append(step_inp["id"])
    out = {"run": run_file,
           "id": name,
           "in": sinputs,
           "out": [{"id": workflow.get_base_id(output["id"])} for output in outputs]}
    if _is_scatter_parallel(parallel):
        assert scatter_inputs, "Did not find items to scatter on: %s" % name
        out.update({"scatterMethod": "dotproduct",
                    "scatter": scatter_inputs})
    return out

def prep_cwl(samples, workflow_fn, out_dir, out_file, integrations=None):
    """Output a CWL description with sub-workflows and steps.
    """
    step_dir = utils.safe_makedir(os.path.join(out_dir, "steps"))
    variables, keyvals = _flatten_samples(samples, out_file, integrations)
    file_estimates = _calc_input_estimates(keyvals, integrations)
    out = _cwl_workflow_template(variables)
    parent_wfs = []
    steps, wfoutputs = workflow_fn(samples)
    used_inputs = set([])
    for cur in workflow.generate(variables, steps, wfoutputs):
        if cur[0] == "step":
            _, name, parallel, inputs, outputs, image, programs, disk, cores = cur
            step_file = _write_tool(step_dir, name, inputs, outputs, parallel, image, programs,
                                    file_estimates, disk, cores, samples)
            out["steps"].append(_step_template(name, step_file, inputs, outputs, parallel))
            used_inputs |= set(x["id"] for x in inputs)
        elif cur[0] == "upload":
            for output in cur[1]:
                wf_output = copy.deepcopy(output)
                if "outputSource" not in wf_output:
                    wf_output["outputSource"] = wf_output.pop("source")
                wf_output = _clean_record(wf_output)
                out["outputs"].append(wf_output)
        elif cur[0] == "wf_start":
            parent_wfs.append(out)
            out = _cwl_workflow_template(cur[1])
        elif cur[0] == "wf_finish":
            _, name, parallel, inputs, outputs, scatter = cur
            wf_out_file = "wf-%s.cwl" % name
            with open(os.path.join(out_dir, wf_out_file), "w") as out_handle:
                yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
            out = parent_wfs.pop(-1)
            out["steps"].append(_step_template(name, wf_out_file, inputs, outputs, parallel, scatter))
            used_inputs |= set(x["id"] for x in inputs)
        else:
            raise ValueError("Unexpected workflow value %s" % str(cur))

    with open(out_file, "w") as out_handle:
        out["inputs"] = [x for x in out["inputs"] if x["id"] in used_inputs]
        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    sample_json = "%s-samples.json" % utils.splitext_plus(out_file)[0]
    out_clean = _clean_final_outputs(copy.deepcopy({k: v for k, v in keyvals.items() if k in used_inputs}),
                                     integrations)
    with open(sample_json, "w") as out_handle:
        json.dump(out_clean, out_handle, sort_keys=True, indent=4, separators=(',', ': '))
    return out_file, sample_json

def _flatten_samples(samples, base_file, integrations=None):
    """Create a flattened JSON representation of data from the bcbio world map.
    """
    flat_data = []
    for data in samples:
        data["reference"] = _indexes_to_secondary_files(data["reference"], data["genome_build"])
        cur_flat = {}
        for key_path in [["analysis"], ["description"], ["rgnames"], ["config", "algorithm"],
                         ["metadata"], ["genome_build"], ["resources"],
                         ["files"], ["reference"], ["genome_resources"], ["vrn_file"]]:
            cur_key = "__".join(key_path)
            for flat_key, flat_val in _to_cwldata(cur_key, tz.get_in(key_path, data)):
                cur_flat[flat_key] = flat_val
        flat_data.append(cur_flat)
    out = {}
    for key in sorted(list(set(reduce(operator.add, [d.keys() for d in flat_data])))):
        # Periods in keys cause issues with WDL and some CWL implementations
        clean_key = key.replace(".", "_")
        out[clean_key] = []
        for cur_flat in flat_data:
            out[clean_key].append(cur_flat.get(key))
    # special case for back-compatibility with fasta specifications -- yuck
    if "reference__fasta__base" not in out and "reference__fasta" in out:
        out["reference__fasta__base"] = out["reference__fasta"]
        del out["reference__fasta"]
    return _samplejson_to_inputs(out), out

def _indexes_to_secondary_files(gresources, genome_build):
    """Convert a list of genome indexes into a single file plus secondary files.

    This ensures that all indices are staged together in a single directory.
    """
    out = {}
    for refname, val in gresources.items():
        if isinstance(val, dict) and "indexes" in val:
            # list of indexes -- aligners
            if len(val.keys()) == 1:
                indexes = val["indexes"]
                if len(indexes) == 0:
                    raise ValueError("Did not find indexes for %s: %s" % (refname, val))
                elif len(indexes) == 1:
                    val = {"indexes": indexes[0]}
                else:
                    val = {"indexes": {"base": indexes[0], "indexes": indexes[1:]}}
            # directory plus indexes -- snpEff
            elif "base" in val and os.path.isdir(val["base"]) and len(val["indexes"]) > 0:
                indexes = val["indexes"]
                val = {"base": indexes[0], "indexes": indexes[1:]}
        elif isinstance(val, dict) and genome_build in val:
            val = _indexes_to_secondary_files(val, genome_build)
        out[refname] = val
    return out

def _add_suppl_info(inp, val):
    """Add supplementary information to inputs from file values.
    """
    inp["type"] = _get_avro_type(val)
    secondary = _get_secondary_files(val)
    if secondary:
        inp["secondaryFiles"] = secondary
    return inp

def _get_secondary_files(val):
    """Retrieve associated secondary files.

    Normalizes input values into definitions of available secondary files.
    """
    out = []
    if isinstance(val, (tuple, list)):
        for x in val:
            for s in _get_secondary_files(x):
                if s and s not in out:
                    out.append(s)
    elif isinstance(val, dict) and (val.get("class") == "File" or "File" in val.get("class")):
        if "secondaryFiles" in val:
            for sf in [x["path"] for x in val["secondaryFiles"]]:
                rext = _get_relative_ext(val["path"], sf)
                if rext and rext not in out:
                    out.append(rext)
    return out

def _get_relative_ext(of, sf):
    """Retrieve relative extension given the original and secondary files.
    """
    def half_finished_trim(orig, prefix):
        return (os.path.basename(prefix).count(".") > 0 and
                os.path.basename(orig).count(".") == os.path.basename(prefix).count("."))
    prefix = os.path.commonprefix([sf, of])
    while prefix.endswith(".") or (half_finished_trim(sf, prefix) and half_finished_trim(of, prefix)):
        prefix = prefix[:-1]
    exts_to_remove = of.replace(prefix, "")
    ext_to_add = sf.replace(prefix, "")
    # Return extensions relative to original
    if not exts_to_remove or exts_to_remove.startswith("."):
        return "^" * exts_to_remove.count(".") + ext_to_add
    else:
        raise ValueError("No cross platform way to reference complex extension: %s %s" % (sf, of))

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
            elif isinstance(ctype, (list, tuple)):
                for x in ctype:
                    if x not in types:
                        types.append(x)
            elif ctype not in types:
                types.append(ctype)
        # handle empty types, allow null or a string "null" sentinel
        if len(types) == 0:
            types = ["null", "string"]
        # collapse arrays for multiple types
        if len(types) > 1 and all(isinstance(t, dict) and t["type"] == "array" for t in types):
            types = [{"type": "array", "items": [t["items"] for t in types]}]
        return {"type": "array", "items": (types[0] if len(types) == 1 else types)}
    elif val is None:
        return ["null", "string"]
    # encode booleans as string True/False and unencode on other side
    elif isinstance(val, bool) or isinstance(val, basestring) and val.lower() in ["true", "false", "none"]:
        return ["string", "null", "boolean"]
    elif isinstance(val, int):
        return "long"
    elif isinstance(val, float):
        return "double"
    else:
        return "string"

def _samplejson_to_inputs(svals):
    """Convert sample output into inputs for CWL configuration files, with types.
    """
    out = []
    for key, val in svals.items():
        out.append(_add_suppl_info({"id": "%s" % key}, val))
    return out

def _to_cwldata(key, val):
    """Convert nested dictionary into CWL data, flatening and marking up files.

    Moves file objects to the top level, enabling insertion in CWL inputs/outputs.
    """
    out = []
    if isinstance(val, dict):
        if len(val) == 2 and "base" in val and "indexes" in val:
            if len(val["indexes"]) > 0 and val["base"] == val["indexes"][0]:
                out.append(("%s__indexes" % key, _item_to_cwldata(val["base"])))
            else:
                out.append((key, _to_cwlfile_with_indexes(val)))
        # Dump shared nested keys like resources as a JSON string
        elif key in workflow.ALWAYS_AVAILABLE:
            out.append((key, _item_to_cwldata(json.dumps(val))))
        else:
            remain_val = {}
            for nkey, nval in val.items():
                cur_nkey = "%s__%s" % (key, nkey)
                cwl_nval = _item_to_cwldata(nval)
                if isinstance(cwl_nval, dict):
                    out.extend(_to_cwldata(cur_nkey, nval))
                elif key in workflow.ALWAYS_AVAILABLE:
                    remain_val[nkey] = nval
                else:
                    out.append((cur_nkey, cwl_nval))
            if remain_val:
                out.append((key, json.dumps(remain_val, sort_keys=True, separators=(',', ':'))))
    else:
        out.append((key, _item_to_cwldata(val)))
    return out

def _to_cwlfile_with_indexes(val):
    """Convert reads with ready to go indexes into the right CWL object.

    Identifies the top level directory and creates a tarball, avoiding
    trying to handle complex secondary setups which are not cross platform.

    Skips doing this for reference files, which take up too much time and
    space to unpack multiple times.
    """
    if val["base"].endswith(".fa") and any([x.endswith(".fa.fai") for x in val["indexes"]]):
        return _item_to_cwldata(val["base"])
    else:
        dirname = os.path.dirname(val["base"])
        assert all([x.startswith(dirname) for x in val["indexes"]])
        return {"class": "File", "path": _directory_tarball(dirname)}

def _item_to_cwldata(x):
    """"Markup an item with CWL specific metadata.
    """
    if isinstance(x, (list, tuple)):
        return [_item_to_cwldata(subx) for subx in x]
    elif (x and isinstance(x, basestring) and
          (((os.path.isfile(x) or os.path.isdir(x)) and os.path.exists(x)) or
           objectstore.is_remote(x))):
        if os.path.isfile(x) or objectstore.is_remote(x):
            out = {"class": "File", "path": x}
            if x.endswith(".bam"):
                out["secondaryFiles"] = [{"class": "File", "path": x + ".bai"}]
            elif x.endswith((".vcf.gz", ".bed.gz")):
                out["secondaryFiles"] = [{"class": "File", "path": x + ".tbi"}]
            elif x.endswith(".fa"):
                secondary = [x + ".fai", os.path.splitext(x)[0] + ".dict"]
                secondary = [y for y in secondary if os.path.exists(y) or objectstore.is_remote(x)]
                if secondary:
                    out["secondaryFiles"] = [{"class": "File", "path": y} for y in secondary]
            elif x.endswith(".fa.gz"):
                secondary = [x + ".fai", x + ".gzi", x.replace(".fa.gz", "") + ".dict"]
                secondary = [y for y in secondary if os.path.exists(y) or objectstore.is_remote(x)]
                if secondary:
                    out["secondaryFiles"] = [{"class": "File", "path": y} for y in secondary]
            elif x.endswith(".fq.gz") or x.endswith(".fastq.gz"):
                secondary = [x + ".gbi"]
                secondary = [y for y in secondary if os.path.exists(y) or objectstore.is_remote(x)]
                if secondary:
                    out["secondaryFiles"] = [{"class": "File", "path": y} for y in secondary]
        else:
            out = {"class": "File", "path": _directory_tarball(x)}
        return out
    elif isinstance(x, bool):
        return str(x)
    else:
        return x

def _directory_tarball(dirname):
    """Create a tarball of a complex directory, avoiding complex secondaryFiles.

    Complex secondary files do not work on multiple platforms and are not portable
    to WDL, so for now we create a tarball that workers will unpack.
    """
    assert os.path.isdir(dirname)
    base_dir, tarball_dir = os.path.split(dirname)
    while base_dir and not os.path.exists(os.path.join(base_dir, "seq")):
        base_dir, extra_tarball = os.path.split(base_dir)
        tarball_dir = os.path.join(extra_tarball, tarball_dir)
    tarball = os.path.join(base_dir, "%s-wf.tar.gz" % (tarball_dir.replace(os.path.sep, "--")))
    if not utils.file_exists(tarball):
        with utils.chdir(base_dir):
            with tarfile.open(tarball, "w:gz") as tar:
                tar.add(tarball_dir)
    return tarball

def _clean_final_outputs(keyvals, integrations=None):
    def clean_path(integrations, x):
        retriever = _get_retriever(x, integrations)
        if retriever:
            return retriever.clean_file(x)
        else:
            return x
    def null_to_string(x):
        """Convert None values into the string 'null'

        Required for platforms like SevenBridges without null support from inputs.
        """
        return "null" if x is None else x
    keyvals = _adjust_items(keyvals, null_to_string)
    keyvals = _adjust_files(keyvals, functools.partial(clean_path, integrations))
    return keyvals

def _adjust_items(xs, adjust_fn):
    if isinstance(xs, (list, tuple)):
        return [_adjust_items(x, adjust_fn) for x in xs]
    elif isinstance(xs, dict):
        out = {}
        for k, v in xs.items():
            out[k] = _adjust_items(v, adjust_fn)
        return out
    else:
        return adjust_fn(xs)

def _adjust_files(xs, adjust_fn):
    """Walk over key/value, tuples applying adjust_fn to files.
    """
    if isinstance(xs, dict):
        if "path" in xs:
            out = {}
            out["path"] = adjust_fn(xs["path"])
            for k, vs in xs.items():
                if k != "path":
                    out[k] = _adjust_files(vs, adjust_fn)
            return out
        else:
            out = {}
            for k, vs in xs.items():
                out[k] = _adjust_files(vs, adjust_fn)
            return out
    elif isinstance(xs, (list, tuple)):
        return [_adjust_files(x, adjust_fn) for x in xs]
    else:
        return xs

def _calc_input_estimates(keyvals, integrations=None):
    """Calculate estimations of input file sizes for disk usage approximation.

    These are current dominated by fastq/BAM sizes, so estimate based on that.
    """
    out = {}
    for key, val in keyvals.items():
        size = _calc_file_size(val, 0, integrations)
        if size:
            out[key] = size
    return out

def _calc_file_size(val, depth, integrations):
    if isinstance(val, (list, tuple)):
        sizes = [_calc_file_size(x, depth + 1, integrations) for x in val]
        sizes = [x for x in sizes if x]
        if sizes:
            # Top level, biggest item, otherwise all files together
            return max(sizes) if depth == 0 else sum(sizes)
    elif isinstance(val, dict) and "path" in val:
        return _get_file_size(val["path"], integrations)
    return None

def _get_retriever(path, integrations):
    if path.startswith(tuple(INTEGRATION_MAP.keys())):
        return integrations.get(INTEGRATION_MAP[path.split(":")[0] + ":"])

def _get_file_size(path, integrations):
    """Return file size in megabytes, including querying remote integrations
    """
    retriever = _get_retriever(path, integrations)
    if retriever:
        return retriever.file_size(path)
    elif os.path.exists(path):
        return os.path.getsize(path) / (1024.0 * 1024.0)
