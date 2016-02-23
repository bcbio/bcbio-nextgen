"""Run distributed functions provided a name and json/YAML file with arguments.

Enables command line access and alternative interfaces to run specific
functionality within bcbio-nextgen.
"""
import json
import os
import pprint

import toolz as tz
import yaml

from bcbio import log, utils
from bcbio.distributed import multitasks
from bcbio.pipeline import config_utils, run_info

def process(args):
    """Run the function in args.name given arguments in args.argfile.
    """
    try:
        fn = getattr(multitasks, args.name)
    except AttributeError:
        raise AttributeError("Did not find exposed function in bcbio.distributed.multitasks named '%s'" % args.name)
    if args.moreargs or args.raw:
        fnargs = [args.argfile] + args.moreargs
        work_dir = None
        argfile = None
    else:
        with open(args.argfile) as in_handle:
            fnargs = yaml.safe_load(in_handle)
        work_dir = os.path.dirname(args.argfile)
        fnargs = config_utils.merge_resources(fnargs)
        argfile = args.outfile if args.outfile else "%s-out%s" % os.path.splitext(args.argfile)
    if not work_dir:
        work_dir = os.getcwd()
    if len(fnargs) > 0 and fnargs[0] == "cwl":
        fnargs, parallel = _world_from_cwl(fnargs[1:], work_dir)
        argfile = os.path.join(work_dir, "cwl-%s-world.json" % args.name)
    with utils.chdir(work_dir):
        log.setup_local_logging(parallel={"wrapper": "runfn"})
        out = fn(fnargs)
    if argfile:
        with open(argfile, "w") as out_handle:
            if argfile.endswith(".json"):
                if parallel in ["single-split", "multi-combined", "batch-split"]:
                    json.dump([utils.to_single_data(xs) for xs in out],
                              out_handle)
                elif parallel in ["multi-batch"]:
                    json.dump([_collapse_to_cwl_record(xs, work_dir) for xs in out], out_handle)
                else:
                    json.dump(utils.to_single_data(utils.to_single_data(out)), out_handle)
            else:
                yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)

def _add_resources(data, runtime):
    if "config" in data:
        data["config"]["resources"] = {"default": {"cores": runtime["cores"],
                                                    "memory": "%sM" % int(float(runtime["ram"]) / runtime["cores"])}}
        data["config"]["algorithm"]["num_cores"] = runtime["cores"]
    return data

def _world_from_cwl(fnargs, work_dir):
    """Reconstitute a bcbio world data object from flattened CWL-compatible inputs.

    Converts the flat CWL representation into a nested bcbio world dictionary.

    Handles single sample inputs (returning a single world object) and multi-sample
    runs (returning a list of individual samples to get processed together).
    """
    parallel = None
    runtime = {}
    out = []
    data = {}
    passed_keys = []
    for fnarg in fnargs:
        key, val = fnarg.split("=")
        # extra values pulling in nested indexes
        if key == "ignore":
            continue
        if key == "sentinel-parallel":
            parallel = val
            continue
        if key == "sentinel-runtime":
            runtime = json.loads(val)
            continue
        # starting a new record -- duplicated key
        if key in passed_keys:
            data["dirs"] = {"work": work_dir}
            data["cwl_keys"] = passed_keys
            data = _add_resources(data, runtime)
            data = run_info.normalize_world(data)
            out.append(data)
            data = {}
            passed_keys = []
        passed_keys.append(key)
        key = key.split("__")
        if val.startswith(("{", "[")):
            val = json.loads(val)
        elif val.find(";;") >= 0:
            val = val.split(";;")
        data = _update_nested(key, val, data)
    if data:
        data["dirs"] = {"work": work_dir}
        data["cwl_keys"] = passed_keys
        data = _add_resources(data, runtime)
        data = run_info.normalize_world(data)
        out.append(data)
    if parallel in ["single-parallel", "single-merge", "multi-parallel", "multi-combined", "multi-batch",
                    "batch-split", "batch-parallel", "batch-merge", "batch-single"]:
        out = [out]
    else:
        assert len(out) == 1, "%s\n%s" % (pprint.pformat(out), pprint.pformat(fnargs))
    return out, parallel

def _collapse_to_cwl_record(samples, work_dir):
    """Convert nested samples from batches into a CWL record, based on input keys.
    """
    all_keys = sorted(list(set().union(*[d["cwl_keys"] for d in samples])), key=lambda x: (-len(x), tuple(x)))
    out = {}
    for key in all_keys:
        key_parts = key.split("__")
        vals = []
        cur = []
        for d in samples:
            val = tz.get_in(key_parts, d)
            if isinstance(val, basestring):
                if os.path.exists(val):
                    val = {"class": "File", "path": val}
                    secondary = []
                    for idx in [".bai", ".tbi", ".gbi", ".fai"]:
                        idx_file = val["path"] + idx
                        if os.path.exists(idx_file):
                            secondary.append({"class": "File", "path": idx_file})
                    for idx in [".dict"]:
                        idx_file = os.path.splitext(val["path"])[0] + idx
                        if os.path.exists(idx_file):
                            secondary.append({"class": "File", "path": idx_file})
                    if secondary:
                        val["secondaryFiles"] = secondary
            elif isinstance(val, dict):
                val = json.dumps(val)
            vals.append(val)
            # Remove nested keys to avoid specifying multiple times
            cur.append(_dissoc_in(d, key_parts) if len(key_parts) > 1 else d)
        samples = cur
        out[key] = vals
    return out

def _dissoc_in(d, key_parts):
    if len(key_parts) > 1:
        d[key_parts[0]] = _dissoc_in(d[key_parts[0]], key_parts[1:])
    else:
        del d[key_parts[0]]
    return d

def _update_nested(key, val, data):
    """Update the data object, avoiding over-writing with nested dictionaries.
    """
    if isinstance(val, dict):
        for sub_key, sub_val in val.items():
            data = _update_nested(key + [sub_key], sub_val, data)
    else:
        if tz.get_in(key, data) is not None:
            raise ValueError("Duplicated key %s" % key)
        data = tz.update_in(data, key, lambda x: val)
    return data

def add_subparser(subparsers):
    parser = subparsers.add_parser("runfn", help=("Run a specific bcbio-nextgen function."
                                                  "Intended for distributed use."))
    parser.add_argument("name", help="Name of the function to run")
    parser.add_argument("argfile", help="JSON file with arguments to the function")
    parser.add_argument("moreargs", nargs="*", help="Additional arguments to pass in the case of raw input")
    parser.add_argument("--raw", action="store_true", default=False,
                        help="Treat the inputs as raw file arguments to the function, instead of parsing them.")
    parser.add_argument("-o", "--outfile", help="Output file to write, defaults to inputfile-out.yaml")
