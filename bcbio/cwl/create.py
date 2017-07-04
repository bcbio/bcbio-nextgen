"""Create Common Workflow Language (CWL) runnable files and tools from a world object.
"""
import copy
import functools
import json
import operator
import os

import toolz as tz
import yaml

from bcbio import utils
from bcbio.cwl import defs, workflow
from bcbio.distributed import objectstore, resources

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
                             {"class": "InlineJavascriptRequirement"},  # for secondary Files
                             {"class": "ScatterFeatureRequirement"},
                             {"class": "SubworkflowFeatureRequirement"}],
            "inputs": ready_inputs,
            "outputs": [],
            "steps": []}

def _write_tool(step_dir, name, inputs, outputs, parallel, image, programs,
                file_estimates, disk, step_cores, samples):
    out_file = os.path.join(step_dir, "%s.cwl" % name)
    resource_cores, mem_gb_per_core = resources.cpu_and_memory((programs or []) + ["default"], samples)
    cores = step_cores if step_cores else resource_cores
    mem_mb_total = int(mem_gb_per_core * cores * 1024)
    bcbio_docker_disk = 1 * 1024  # Minimum requirements for bcbio Docker image
    cwl_res = {"class": "ResourceRequirement",
               "coresMin": cores, "ramMin": mem_mb_total, "outdirMin": bcbio_docker_disk}
    docker_image = "bcbio/bcbio" if image == "bcbio" else "quay.io/bcbio/%s" % image
    docker = {"class": "DockerRequirement", "dockerPull": docker_image, "dockerImageId": docker_image}
    if file_estimates and disk:
        total_estimate = 0
        for key, multiplier in disk.items():
            if key in file_estimates:
                total_estimate += int(multiplier * file_estimates[key])
        if total_estimate:
            cwl_res["tmpdirMin"] = total_estimate
            cwl_res["outdirMin"] += total_estimate
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

def _get_sentinel_val(v):
    """Retrieve expected sentinel value for an output, expanding records.
    """
    out = workflow.get_base_id(v["id"])
    if workflow.is_cwl_record(v):
        def _get_fields(d):
            if isinstance(d, dict):
                if "fields" in d:
                    return d["fields"]
                else:
                    for v in d.values():
                        fields = _get_fields(v)
                        if fields:
                            return fields
        out += ":%s" % ";".join([x["name"] for x in _get_fields(v)])
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
    sample_json, variables, keyvals = _flatten_samples(samples, out_file, integrations)
    file_estimates = _calc_input_estimates(keyvals, integrations)
    out = _cwl_workflow_template(variables)
    parent_wfs = []
    steps, wfoutputs = workflow_fn(samples)
    for cur in workflow.generate(variables, steps, wfoutputs):
        if cur[0] == "step":
            _, name, parallel, inputs, outputs, image, programs, disk, cores = cur
            step_file = _write_tool(step_dir, name, inputs, outputs, parallel, image, programs,
                                    file_estimates, disk, cores, samples)
            out["steps"].append(_step_template(name, step_file, inputs, outputs, parallel))
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
        else:
            raise ValueError("Unexpected workflow value %s" % str(cur))

    with open(out_file, "w") as out_handle:
        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out_file, sample_json

def _flatten_samples(samples, base_file, integrations=None):
    """Create a flattened JSON representation of data from the bcbio world map.
    """
    out_file = "%s-samples.json" % utils.splitext_plus(base_file)[0]
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
        out[key] = []
        for cur_flat in flat_data:
            out[key].append(cur_flat.get(key))
    # special case for back-compatibility with fasta specifications -- yuck
    if "reference__fasta__base" not in out and "reference__fasta" in out:
        out["reference__fasta__base"] = out["reference__fasta"]
        del out["reference__fasta"]
    out_clean = _clean_final_outputs(copy.deepcopy(out), integrations)
    with open(out_file, "w") as out_handle:
        json.dump(out_clean, out_handle, sort_keys=True, indent=4, separators=(',', ': '))
    return out_file, _samplejson_to_inputs(out), out

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
                if len(indexes) == 1:
                    val = {"indexes": indexes[0]}
                else:
                    val = {"indexes": {"base": indexes[0], "indexes": indexes[1:]}}
            # directory plus indexes -- snpEff
            elif "base" in val and os.path.isdir(val["base"]) and len(val["indexes"]) > 0:
                indexes = val["indexes"]
                val = {"base": indexes[0]}
                if len(indexes) > 1:
                    val["indexes"] = indexes[1:]
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
    # Return entire file relative to original
    # no way to cleanly reference dirname without using InlineJavascriptRequirement
    elif prefix.endswith("/"):
        return '$(self.location.substr(0, self.location.lastIndexOf("/")))/%s' % ext_to_add

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
    elif isinstance(val, bool):
        return "string"
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
            out.append((key, _to_cwlfile_with_indexes(val)))
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
    """
    return {"class": "File", "path": val["base"],
            "secondaryFiles": [{"class": "File", "path": x} for x in val["indexes"]]}

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
        else:
            # aligner and database indices where we list the entire directory as secondary files
            dir_targets = ("mainIndex", ".alt", ".amb", ".ann", ".bwt", ".pac", ".sa", ".ebwt", ".bt2",
                           "Genome", "GenomeIndex", "GenomeIndexHash", "OverflowTable")
            assert os.path.isdir(x)
            base_name = None
            fnames = sorted(os.listdir(x))
            for fname in fnames:
                if fname.endswith(dir_targets):
                    base_name = fname
                    break
            if base_name:
                fnames.pop(fnames.index(base_name))
                base_name = os.path.join(x, base_name)
                fnames = [os.path.join(x, y) for y in fnames]
                out = {"class": "File", "path": base_name,
                       "secondaryFiles": [{"class": "File", "path": f} for f in fnames]}
            # skip directories we're not currently using in CWL recipes
            else:
                out = None
        return out
    elif isinstance(x, bool):
        return str(x)
    else:
        return x

def _clean_final_outputs(keyvals, integrations=None):
    def clean_path(integrations, x):
        retriever = _get_retriever(x, integrations)
        if retriever:
            return retriever.clean_file(x)
        else:
            return x
    return _adjust_files(keyvals, functools.partial(clean_path, integrations))

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
    input_sizes = []
    for sample_files in keyvals["files"]:
        sample_inputs = 0
        for cwl_file in sample_files:
            file_size = _get_file_size(cwl_file["path"], integrations)
            if file_size:
                sample_inputs += file_size
        if sample_inputs > 0:
            input_sizes.append(sample_inputs)
    if len(input_sizes) > 0:
        return {"files": max(input_sizes)}

def _get_retriever(path, integrations):
    integration_map = {"keep:": "arvados", "s3:": "s3", "sbg:": "sbgenomics"}
    if path.startswith(tuple(integration_map.keys())):
        return integrations.get(integration_map[path.split(":")[0] + ":"])

def _get_file_size(path, integrations):
    """Return file size in megabytes, including querying remote integrations
    """
    retriever = _get_retriever(path, integrations)
    if retriever:
        return retriever.file_size(path)
    elif os.path.exists(path):
        return os.path.getsize(path) / (1024.0 * 1024.0)
