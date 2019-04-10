"""Create Common Workflow Language (CWL) runnable files and tools from a world object.
"""
from __future__ import print_function
import collections
import copy
import dateutil
import functools
import json
import math
import operator
import os
import tarfile

import requests
import six
import toolz as tz
import yaml

from bcbio import utils
from bcbio.cwl import defs, workflow
from bcbio.distributed import objectstore, resources
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import alignment
from functools import reduce

INTEGRATION_MAP = {"keep:": "arvados", "s3:": "s3", "sbg:": "sbgenomics",
                   "dx:": "dnanexus", "gs:": "gs"}

def from_world(world, run_info_file, integrations=None, add_container_tag=None):
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
    prep_cwl(samples, workflow_fn, out_dir, out_file, integrations, add_container_tag=add_container_tag)

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

def _get_disk_estimates(name, parallel, inputs, file_estimates, samples, disk,
                        cur_remotes, no_files):
    """Retrieve disk usage estimates as CWL ResourceRequirement and hint.

    Disk specification for temporary files and outputs.

    Also optionally includes disk input estimates as a custom hint for
    platforms which need to stage these and don't pre-estimate these when
    allocating machine sizes.
    """
    tmp_disk, out_disk, in_disk = 0, 0, 0
    if file_estimates:
        if disk:
            for key, multiplier in disk.items():
                if key in file_estimates:
                    out_disk += int(multiplier * file_estimates[key])
        for inp in inputs:
            scale = 2.0 if inp.get("type") == "array" else 1.0
            # Allocating all samples, could remove for `to_rec` when we ensure we
            # don't have to stage. Currently dnanexus stages everything so need to consider
            if parallel in ["multi-combined", "multi-batch"] and "dnanexus" in cur_remotes:
                scale *= (len(samples))
            if workflow.is_cwl_record(inp):
                for f in _get_record_fields(inp):
                    if f["name"] in file_estimates:
                        in_disk += file_estimates[f["name"]] * scale
            elif inp["id"] in file_estimates:
                in_disk += file_estimates[inp["id"]] * scale
        # Round total estimates to integer, assign extra half to temp space
        # It's not entirely clear how different runners interpret this
        tmp_disk = int(math.ceil(out_disk * 0.5))
        out_disk = int(math.ceil(out_disk))

    bcbio_docker_disk = (10 if cur_remotes else 1) * 1024  # Minimum requirements for bcbio Docker image
    disk_hint = {"outdirMin": bcbio_docker_disk + out_disk, "tmpdirMin": tmp_disk}
    # Skip input disk for steps which require only transformation (and thus no staging)
    if no_files:
        in_disk = 0
    # Avoid accidentally flagging as no staging if we don't know sizes of expected inputs
    elif in_disk == 0:
        in_disk = 1
    input_hint = {"class": "dx:InputResourceRequirement", "indirMin": int(math.ceil(in_disk))}
    return disk_hint, input_hint

def _add_current_quay_tag(repo, container_tags):
    """Lookup the current quay tag for the repository, adding to repo string.

    Enables generation of CWL explicitly tied to revisions.
    """
    if ':' in repo:
        return repo, container_tags
    try:
        latest_tag = container_tags[repo]
    except KeyError:
        repo_id = repo[repo.find('/') + 1:]
        tags = requests.request("GET", "https://quay.io/api/v1/repository/" + repo_id).json()["tags"]
        latest_tag = None
        latest_modified = None
        for tag, info in tags.items():
            if latest_tag:
                if (dateutil.parser.parse(info['last_modified']) > dateutil.parser.parse(latest_modified)
                      and tag != 'latest'):
                    latest_modified = info['last_modified']
                    latest_tag = tag
            else:
                latest_modified = info['last_modified']
                latest_tag = tag
        container_tags[repo] = str(latest_tag)
    latest_pull = repo + ':' + str(latest_tag)
    return latest_pull, container_tags

def _write_tool(step_dir, name, inputs, outputs, parallel, image, programs,
                file_estimates, disk, step_cores, samples, cur_remotes, no_files,
                container_tags=None):
    out_file = os.path.join(step_dir, "%s.cwl" % name)
    resource_cores, mem_gb_per_core = resources.cpu_and_memory((programs or []) + ["default"], samples)
    cores = min([step_cores, resource_cores]) if step_cores else resource_cores
    mem_mb_total = int(mem_gb_per_core * cores * 1024)
    cwl_res = {"class": "ResourceRequirement", "coresMin": cores, "ramMin": mem_mb_total}
    disk_hint, input_hint = _get_disk_estimates(name, parallel, inputs, file_estimates, samples, disk,
                                                cur_remotes, no_files)
    cwl_res.update(disk_hint)
    docker_image = "bcbio/bcbio" if image == "bcbio" else "quay.io/bcbio/%s" % image
    if container_tags is not None:
        docker_image, container_tags = _add_current_quay_tag(docker_image, container_tags)
    docker = {"class": "DockerRequirement", "dockerPull": docker_image, "dockerImageId": docker_image}
    out = {"class": "CommandLineTool",
           "cwlVersion": "v1.0",
           "baseCommand": ["bcbio_nextgen.py", "runfn", name, "cwl"],
           "requirements": [],
           "hints": [docker, cwl_res, input_hint],
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
        # GATK requires networking for setting up log4j logging, use arvados extension
        if any(p.startswith(("gatk", "sentieon")) for p in programs):
            out["hints"] += [{"class": "arv:APIRequirement"}]
    # Multi-process methods that read heavily from BAM files need extra keep cache for Arvados
    if name in ["pipeline_summary", "variantcall_batch_region", "detect_sv"]:
        out["hints"] += [{"class": "arv:RuntimeConstraints", "keep_cache": 4096}]
    def add_to_namespaces(k, v, out):
        if "$namespaces" not in out:
            out["$namespaces"] = {}
        out["$namespaces"][k] = v
        return out
    if any(h.get("class", "").startswith("arv:") for h in out["hints"]):
        out = add_to_namespaces("arv", "http://arvados.org/cwl#", out)
    if any(h.get("class", "").startswith("dx") for h in out["hints"]):
        out = add_to_namespaces("dx", "https://www.dnanexus.com/cwl#", out)
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
                                                          for v in inputs]),
                         "run_number=0"]
    out = _add_inputs_to_tool(inputs, out, parallel, use_commandline_args)
    out = _add_outputs_to_tool(outputs, out)
    _tool_to_file(out, out_file)
    return os.path.join("steps", os.path.basename(out_file))

def _write_expressiontool(step_dir, name, inputs, outputs, expression, parallel):
    """Create an ExpressionTool output for the given inputs
    """
    out_file = os.path.join(step_dir, "%s.cwl" % name)
    out = {"class": "ExpressionTool",
           "cwlVersion": "v1.0",
           "requirements": [{"class": "InlineJavascriptRequirement"}],
           "inputs": [],
           "outputs": [],
           "expression": expression}
    out = _add_inputs_to_tool(inputs, out, parallel)
    out = _add_outputs_to_tool(outputs, out)
    _tool_to_file(out, out_file)
    return os.path.join("steps", os.path.basename(out_file))

def _add_outputs_to_tool(outputs, tool):
    for outp in outputs:
        outp_tool = copy.deepcopy(outp)
        outp_tool = _clean_record(outp_tool)
        outp_tool["id"] = workflow.get_base_id(outp["id"])
        tool["outputs"].append(outp_tool)
    return tool

def _add_inputs_to_tool(inputs, tool, parallel, use_commandline_args=False):
    for i, inp in enumerate(inputs):
        base_id = workflow.get_base_id(inp["id"])
        inp_tool = copy.deepcopy(inp)
        inp_tool["id"] = base_id
        if inp.get("wf_duplicate"):
            inp_tool["id"] += "_toolinput"
        for attr in ["source", "valueFrom", "wf_duplicate"]:
            inp_tool.pop(attr, None)
        # Ensure records and workflow inputs get scattered
        if (_is_scatter_parallel(parallel) and _do_scatter_var(inp, parallel) and
              (workflow.is_cwl_record(inp) or inp["wf_duplicate"])):
            inp_tool = workflow._flatten_nested_input(inp_tool)
        if use_commandline_args:
            inp_binding = {"prefix": "%s=" % base_id,
                           "separate": False, "itemSeparator": ";;", "position": i}
            inp_tool = _place_input_binding(inp_tool, inp_binding, parallel)
        else:
            inp_binding = None
        inp_tool = _place_secondary_files(inp_tool, inp_binding)
        inp_tool = _clean_record(inp_tool)
        tool["inputs"].append(inp_tool)
    return tool

def _tool_to_file(tool, out_file):
    with open(out_file, "w") as out_handle:
        def str_presenter(dumper, data):
            if len(data.splitlines()) > 1:  # check for multiline string
                return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')
            return dumper.represent_scalar('tag:yaml.org,2002:str', data)
        yaml.add_representer(str, str_presenter)
        yaml.dump(tool, out_handle, default_flow_style=False, allow_unicode=False)
    return out_file

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
        return (val == "File" or (isinstance(val, (list, tuple)) and
                                  ("File" in val or any(isinstance(x, dict) and _is_file(val)) for x in val)))
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

def _step_template(name, run_file, inputs, outputs, parallel, step_parallelism, scatter=None):
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
        # An initial parallel scatter and multiple chained parallel sample scatters
        if (parallel == "multi-parallel" and
              (not step_parallelism or
               step_parallelism.get(workflow.get_step_prefix(inp["id"])) == "multi-parallel")):
            scatter_inputs.append(step_inp["id"])
        # scatter on inputs from previous processes that have been arrayed
        elif (_is_scatter_parallel(parallel) and (_do_scatter_var(inp, parallel)
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

def _get_cur_remotes(path):
    """Retrieve remote references defined in the CWL.
    """
    cur_remotes = set([])
    if isinstance(path, (list, tuple)):
        for v in path:
            cur_remotes |= _get_cur_remotes(v)
    elif isinstance(path, dict):
        for v in path.values():
            cur_remotes |= _get_cur_remotes(v)
    elif path and isinstance(path, six.string_types):
        if path.startswith(tuple(INTEGRATION_MAP.keys())):
            cur_remotes.add(INTEGRATION_MAP.get(path.split(":")[0] + ":"))
    return cur_remotes

def prep_cwl(samples, workflow_fn, out_dir, out_file, integrations=None,
             add_container_tag=None):
    """Output a CWL description with sub-workflows and steps.
    """
    if add_container_tag is None:
        container_tags = None
    elif add_container_tag.lower() == "quay_lookup":
        container_tags = {}
    else:
        container_tags = collections.defaultdict(lambda: add_container_tag)
    step_dir = utils.safe_makedir(os.path.join(out_dir, "steps"))
    get_retriever = GetRetriever(integrations, samples)
    variables, keyvals = _flatten_samples(samples, out_file, get_retriever)
    cur_remotes = _get_cur_remotes(keyvals)
    file_estimates = _calc_input_estimates(keyvals, get_retriever)
    out = _cwl_workflow_template(variables)
    parent_wfs = []
    step_parallelism = {}
    steps, wfoutputs = workflow_fn(samples)
    used_inputs = set([])
    for cur in workflow.generate(variables, steps, wfoutputs):
        if cur[0] == "step":
            _, name, parallel, inputs, outputs, image, programs, disk, cores, no_files = cur
            step_file = _write_tool(step_dir, name, inputs, outputs, parallel, image, programs,
                                    file_estimates, disk, cores, samples, cur_remotes, no_files, container_tags)
            out["steps"].append(_step_template(name, step_file, inputs, outputs, parallel, step_parallelism))
            used_inputs |= set(x["id"] for x in inputs)
        elif cur[0] == "expressiontool":
            _, name, inputs, outputs, expression, parallel = cur
            step_file = _write_expressiontool(step_dir, name, inputs, outputs, expression, parallel)
            out["steps"].append(_step_template(name, step_file, inputs, outputs, parallel, step_parallelism))
            used_inputs |= set(x["id"] for x in inputs)
        elif cur[0] == "upload":
            for output in cur[1]:
                wf_output = copy.deepcopy(output)
                if "outputSource" not in wf_output:
                    wf_output["outputSource"] = wf_output.pop("source")
                wf_output = _clean_record(wf_output)
                # Avoid input/output naming clashes
                if wf_output["id"] in used_inputs:
                    wf_output["id"] = "%s_out" % wf_output["id"]
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
            out["steps"].append(_step_template(name, wf_out_file, inputs, outputs, parallel,
                                               step_parallelism, scatter))
            used_inputs |= set(x["id"] for x in inputs)
        else:
            raise ValueError("Unexpected workflow value %s" % str(cur))
        step_parallelism[name] = parallel

    with open(out_file, "w") as out_handle:
        out["inputs"] = [x for x in out["inputs"] if x["id"] in used_inputs]
        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    sample_json = "%s-samples.json" % utils.splitext_plus(out_file)[0]
    out_clean = _clean_final_outputs(copy.deepcopy({k: v for k, v in keyvals.items() if k in used_inputs}),
                                     get_retriever)
    with open(sample_json, "w") as out_handle:
        json.dump(out_clean, out_handle, sort_keys=True, indent=4, separators=(',', ': '))
    return out_file, sample_json

def _flatten_samples(samples, base_file, get_retriever):
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
            for flat_key, flat_val in _to_cwldata(cur_key, tz.get_in(key_path, data), get_retriever):
                cur_flat[flat_key] = flat_val
        flat_data.append(cur_flat)
    out = {}
    for key in sorted(list(set(reduce(operator.add, [list(d.keys()) for d in flat_data])))):
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
                indexes = sorted(val["indexes"])
                if len(indexes) == 0:
                    if refname not in alignment.allow_noindices():
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
    Requires indices to be present in all files, since declared CWL secondary
    files are not optional. So if we have a mix of BAM (bai) and fastq (gbi) we
    ignore the existing indices and will have to regenerate during processing.
    """
    out = []
    if isinstance(val, (tuple, list)):
        s_counts = collections.defaultdict(int)
        for x in val:
            for s in _get_secondary_files(x):
                s_counts[s] += 1
        for s, count in s_counts.items():
            if s and s not in out and count == len([x for x in val if x]):
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
    # Handle remote files
    if of.find(":") > 0:
        of = os.path.basename(of.split(":")[-1])
    if sf.find(":") > 0:
        sf = os.path.basename(sf.split(":")[-1])
    prefix = os.path.commonprefix([sf, of])
    while prefix.endswith(".") or (half_finished_trim(sf, prefix) and half_finished_trim(of, prefix)):
        prefix = prefix[:-1]
    exts_to_remove = of.replace(prefix, "")
    ext_to_add = sf.replace(prefix, "")
    # Return extensions relative to original
    if not exts_to_remove or exts_to_remove.startswith("."):
        return str("^" * exts_to_remove.count(".") + ext_to_add)
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
                    if isinstance(ctype["items"], (list, tuple)):
                        for t in ctype["items"]:
                            if t not in types:
                                types.append(t)
                    else:
                        if ctype not in types:
                            types.append(ctype)
            elif isinstance(ctype, (list, tuple)):
                for x in ctype:
                    if x not in types:
                        types.append(x)
            elif ctype not in types:
                types.append(ctype)
        # handle empty types, allow null
        if len(types) == 0:
            types = ["null"]
            # empty lists
            if isinstance(val, (list, tuple)) and len(val) == 0:
                types.append({"type": "array", "items": ["null"]})
        types = _avoid_duplicate_arrays(types)
        # Avoid empty null only arrays which confuse some runners
        if len(types) == 1 and types[0] == "null":
            types.append("string")
        return {"type": "array", "items": (types[0] if len(types) == 1 else types)}
    elif val is None:
        return ["null"]
    # encode booleans as string True/False and unencode on other side
    elif isinstance(val, bool) or isinstance(val, six.string_types) and val.lower() in ["true", "false", "none"]:
        return ["string", "null", "boolean"]
    elif isinstance(val, int):
        return "long"
    elif isinstance(val, float):
        return "double"
    else:
        return "string"

def _avoid_duplicate_arrays(types):
    """Collapse arrays when we have multiple types.
    """
    arrays = [t for t in types if isinstance(t, dict) and t["type"] == "array"]
    others = [t for t in types if not (isinstance(t, dict) and t["type"] == "array")]
    if arrays:
        items = set([])
        for t in arrays:
            if isinstance(t["items"], (list, tuple)):
                items |= set(t["items"])
            else:
                items.add(t["items"])
        if len(items) == 1:
            items = items.pop()
        else:
            items = sorted(list(items))
        arrays = [{"type": "array", "items": items}]
    return others + arrays

def _samplejson_to_inputs(svals):
    """Convert sample output into inputs for CWL configuration files, with types.
    """
    out = []
    for key, val in svals.items():
        out.append(_add_suppl_info({"id": "%s" % key}, val))
    return out

def _to_cwldata(key, val, get_retriever):
    """Convert nested dictionary into CWL data, flatening and marking up files.

    Moves file objects to the top level, enabling insertion in CWL inputs/outputs.
    """
    out = []
    if isinstance(val, dict):
        if len(val) == 2 and "base" in val and "indexes" in val:
            if len(val["indexes"]) > 0 and val["base"] == val["indexes"][0]:
                out.append(("%s__indexes" % key, _item_to_cwldata(val["base"], get_retriever)))
            else:
                out.append((key, _to_cwlfile_with_indexes(val, get_retriever)))
        # Dump shared nested keys like resources as a JSON string
        elif key in workflow.ALWAYS_AVAILABLE or key in workflow.STRING_DICT:
            out.append((key, _item_to_cwldata(json.dumps(val), get_retriever)))
        elif key in workflow.FLAT_DICT:
            flat = []
            for k, vs in val.items():
                if not isinstance(vs, (list, tuple)):
                    vs = [vs]
                for v in vs:
                    flat.append("%s:%s" % (k, v))
            out.append((key, _item_to_cwldata(flat, get_retriever)))
        else:
            remain_val = {}
            for nkey, nval in val.items():
                cur_nkey = "%s__%s" % (key, nkey)
                cwl_nval = _item_to_cwldata(nval, get_retriever)
                if isinstance(cwl_nval, dict):
                    out.extend(_to_cwldata(cur_nkey, nval, get_retriever))
                elif key in workflow.ALWAYS_AVAILABLE:
                    remain_val[nkey] = nval
                else:
                    out.append((cur_nkey, cwl_nval))
            if remain_val:
                out.append((key, json.dumps(remain_val, sort_keys=True, separators=(',', ':'))))
    else:
        out.append((key, _item_to_cwldata(val, get_retriever)))
    return out

def _remove_remote_prefix(f):
    """Remove any remote references to allow object lookups by file paths.
    """
    return f.split(":")[-1].split("/", 1)[1] if objectstore.is_remote(f) else f

def _index_blacklist(xs):
    blacklist = ["-resources.yaml"]
    return [x for x in xs if not any([x.find(b) >=0 for b in blacklist])]

def _to_cwlfile_with_indexes(val, get_retriever):
    """Convert reads with ready to go indexes into the right CWL object.

    Identifies the top level directory and creates a tarball, avoiding
    trying to handle complex secondary setups which are not cross platform.

    Skips doing this for reference files and standard setups like bwa, which
    take up too much time and space to unpack multiple times.
    """
    val["indexes"] = _index_blacklist(val["indexes"])
    tval = {"base": _remove_remote_prefix(val["base"]),
            "indexes": [_remove_remote_prefix(f) for f in val["indexes"]]}
    # Standard named set of indices, like bwa
    # Do not include snpEff, which we need to isolate inside a nested directory
    # hisat2 indices do also not localize cleanly due to compilicated naming
    cp_dir, cp_base = os.path.split(os.path.commonprefix([tval["base"]] + tval["indexes"]))
    if (cp_base and cp_dir == os.path.dirname(tval["base"]) and
            not ("/snpeff/" in cp_dir or "/hisat2" in cp_dir)):
        return _item_to_cwldata(val["base"], get_retriever, val["indexes"])
    else:
        dirname = os.path.dirname(tval["base"])
        assert all([x.startswith(dirname) for x in tval["indexes"]])
        return {"class": "File", "path": directory_tarball(dirname)}

def _add_secondary_if_exists(secondary, out, get_retriever):
    """Add secondary files only if present locally or remotely.
    """
    secondary = [_file_local_or_remote(y, get_retriever) for y in secondary]
    secondary = [z for z in secondary if z]
    if secondary:
        out["secondaryFiles"] = [{"class": "File", "path": f} for f in secondary]
    return out

def _item_to_cwldata(x, get_retriever, indexes=None):
    """"Markup an item with CWL specific metadata.
    """
    if isinstance(x, (list, tuple)):
        return [_item_to_cwldata(subx, get_retriever) for subx in x]
    elif (x and isinstance(x, six.string_types) and
          (((os.path.isfile(x) or os.path.isdir(x)) and os.path.exists(x)) or
           objectstore.is_remote(x))):
        if _file_local_or_remote(x, get_retriever):
            out = {"class": "File", "path": x}
            if indexes:
                out = _add_secondary_if_exists(indexes, out, get_retriever)
            elif x.endswith(".bam"):
                out = _add_secondary_if_exists([x + ".bai"], out, get_retriever)
            elif x.endswith(".cram"):
                out = _add_secondary_if_exists([x + ".crai"], out, get_retriever)
            elif x.endswith((".vcf.gz", ".bed.gz")):
                out = _add_secondary_if_exists([x + ".tbi"], out, get_retriever)
            elif x.endswith(".fa"):
                out = _add_secondary_if_exists([x + ".fai", os.path.splitext(x)[0] + ".dict"], out, get_retriever)
            elif x.endswith(".fa.gz"):
                out = _add_secondary_if_exists([x + ".fai", x + ".gzi", x.replace(".fa.gz", "") + ".dict"],
                                               out, get_retriever)
            elif x.endswith(".fq.gz") or x.endswith(".fastq.gz"):
                out = _add_secondary_if_exists([x + ".gbi"], out, get_retriever)
            elif x.endswith(".gtf"):
                out = _add_secondary_if_exists([x + ".db"], out, get_retriever)
        else:
            out = {"class": "File", "path": directory_tarball(x)}
        return out
    elif isinstance(x, bool):
        return str(x)
    else:
        return x

def _file_local_or_remote(f, get_retriever):
    """Check for presence of a local or remote file.
    """
    if os.path.exists(f):
        return f
    integration, config = get_retriever.integration_and_config(f)
    if integration:
        return integration.file_exists(f, config)

def directory_tarball(dirname):
    """Create a tarball of a complex directory, avoiding complex secondaryFiles.

    Complex secondary files do not work on multiple platforms and are not portable
    to WDL, so for now we create a tarball that workers will unpack.
    """
    assert os.path.isdir(dirname), dirname
    base_dir, tarball_dir = os.path.split(dirname)
    while not os.path.exists(os.path.join(base_dir, "seq")) and base_dir and base_dir != "/":
        base_dir, extra_tarball = os.path.split(base_dir)
        tarball_dir = os.path.join(extra_tarball, tarball_dir)
    if base_dir == "/" and not os.path.exists(os.path.join(base_dir, "seq")):
        raise ValueError("Did not find relative directory to create tarball for %s" % dirname)
    tarball = os.path.join(base_dir, "%s-wf.tar.gz" % (tarball_dir.replace(os.path.sep, "--")))
    if not utils.file_exists(tarball):
        print("Preparing CWL input tarball: %s" % tarball)
        with file_transaction({}, tarball) as tx_tarball:
            with utils.chdir(base_dir):
                with tarfile.open(tx_tarball, "w:gz") as tar:
                    tar.add(tarball_dir)
    return tarball

def _clean_final_outputs(keyvals, get_retriever):
    def clean_path(get_retriever, x):
        integration, config = get_retriever.integration_and_config(x)
        if integration:
            return integration.clean_file(x, config)
        else:
            return x
    keyvals = _adjust_files(keyvals, functools.partial(clean_path, get_retriever))
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

def _calc_input_estimates(keyvals, get_retriever):
    """Calculate estimations of input file sizes for disk usage approximation.

    These are current dominated by fastq/BAM sizes, so estimate based on that.
    """
    out = {}
    for key, val in keyvals.items():
        size = _calc_file_size(val, 0, get_retriever)
        if size:
            out[key] = size
    return out

def _calc_file_size(val, depth, get_retriever):
    if isinstance(val, (list, tuple)):
        sizes = [_calc_file_size(x, depth + 1, get_retriever) for x in val]
        sizes = [x for x in sizes if x]
        if sizes:
            # Top level, biggest item, otherwise all files together
            return max(sizes) if depth == 0 else sum(sizes)
    elif isinstance(val, dict) and "path" in val:
        return _get_file_size(val["path"], get_retriever)
    return None

class GetRetriever:
    def __init__(self, integrations, samples):
        self._integrations = integrations
        self._samples = samples

    def integration_and_config(self, path):
        """Get a retriever and configuration for the given file path.
        """
        if path.startswith(tuple(INTEGRATION_MAP.keys())):
            key = INTEGRATION_MAP[path.split(":")[0] + ":"]
            integration = self._integrations.get(key)
            config = {}
            for sample in self._samples:
                config = tz.get_in(["config", key], sample)
                if config:
                    break
            return integration, config

        return None, None

def _get_file_size(path, get_retriever):
    """Return file size in megabytes, including querying remote integrations
    """
    integration, config = get_retriever.integration_and_config(path)
    if integration:
        return integration.file_size(path, config)
    elif os.path.exists(path):
        return os.path.getsize(path) / (1024.0 * 1024.0)
