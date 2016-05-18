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
from bcbio.log import logger
from bcbio.distributed import multitasks
from bcbio.pipeline import config_utils, run_info

def process(args):
    """Run the function in args.name given arguments in args.argfile.
    """
    # Set environment to standard to use periods for decimals and avoid localization
    os.environ["LC_ALL"] = "C"
    os.environ["LC"] = "C"
    os.environ["LANG"] = "C"
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
        fnargs = config_utils.merge_resources(fnargs)
        argfile = os.path.join(work_dir, "cwl.output.json")
    else:
        parallel = None
    with utils.chdir(work_dir):
        log.setup_local_logging(parallel={"wrapper": "runfn"})
        try:
            out = fn(fnargs)
        except:
            logger.exception()
            raise
    if argfile:
        try:
            _write_out_argfile(argfile, out, fnargs, parallel, work_dir)
        except:
            logger.exception()
            raise

def _write_out_argfile(argfile, out, fnargs, parallel, work_dir):
    """Write output argfile, preparing a CWL ready JSON or YAML representation of the world.
    """
    with open(argfile, "w") as out_handle:
        if argfile.endswith(".json"):
            if parallel in ["single-split", "multi-combined", "batch-split"]:
                json.dump(_convert_to_cwl_json([utils.to_single_data(xs) for xs in out], fnargs),
                            out_handle, sort_keys=True, indent=4, separators=(', ', ': '))
            elif parallel in ["multi-batch"]:
                json.dump(_combine_cwl_records([_collapse_to_cwl_record(xs, work_dir) for xs in out],
                                               fnargs),
                            out_handle, sort_keys=True, indent=4, separators=(', ', ': '))
            else:
                json.dump(_convert_to_cwl_json(utils.to_single_data(utils.to_single_data(out)), fnargs),
                            out_handle, sort_keys=True, indent=4, separators=(', ', ': '))
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
    output_cwl_keys = None
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
        if key == "sentinel-outputs":
            output_cwl_keys = json.loads(val)
            continue
        # starting a new record -- duplicated key
        if key in passed_keys:
            out.append(_finalize_cwl_in(data, work_dir, passed_keys, output_cwl_keys, runtime))
            data = {}
            passed_keys = []
        passed_keys.append(key)
        key = key.split("__")
        data = _update_nested(key, _convert_value(val), data)
    if data:
        out.append(_finalize_cwl_in(data, work_dir, passed_keys, output_cwl_keys, runtime))
    if parallel in ["single-parallel", "single-merge", "multi-parallel", "multi-combined", "multi-batch",
                    "batch-split", "batch-parallel", "batch-merge", "batch-single"]:
        out = [out]
    else:
        assert len(out) == 1, "%s\n%s" % (pprint.pformat(out), pprint.pformat(fnargs))
    return out, parallel

def _finalize_cwl_in(data, work_dir, passed_keys, output_cwl_keys, runtime):
    """Finalize data object with inputs from CWL.
    """
    data["dirs"] = {"work": work_dir}
    if not tz.get_in(["config", "algorithm"], data):
        if "config" not in data:
            data["config"] = {}
        data["config"]["algorithm"] = {}
    if "rgnames" not in data and "description" in data:
        data["rgnames"] = {"sample": data["description"]}
    data["cwl_keys"] = passed_keys
    data["output_cwl_keys"] = output_cwl_keys
    data = _add_resources(data, runtime)
    data = run_info.normalize_world(data)
    return data

def _convert_value(val):
    """Handle multiple input type values.
    """
    def _is_number(x, op):
        try:
            op(x)
            return True
        except ValueError:
            return False
    if _is_number(val, int):
        return int(val)
    elif _is_number(val, float):
        return float(val)
    elif val.startswith(("{", "[")):
        return json.loads(val)
    elif val.find(";;") >= 0:
        return [_convert_value(v) for v in val.split(";;")]
    elif val.lower() == "true":
        return True
    elif val.lower() == "false":
        return False
    else:
        return val

def _convert_to_cwl_json(data, fnargs):
    """Convert world data object (or list of data objects) into outputs for CWL ingestion.
    """
    out = {}
    for outvar in _get_output_cwl_keys(fnargs):
        keys = outvar.split("__")
        if isinstance(data, dict):
            out[outvar] = _to_cwl(tz.get_in(keys, data))
        else:
            out[outvar] = [_to_cwl(tz.get_in(keys, x)) for x in data]
    return out

def _get_output_cwl_keys(fnargs):
    """Retrieve output_cwl_keys from potentially nested input arguments.
    """
    for items in fnargs:
        if isinstance(items, dict):
            items = [items]
        for d in items:
            if isinstance(d, dict) and d.get("output_cwl_keys"):
                return d["output_cwl_keys"]
    raise ValueError("Did not find output_cwl_keys in %s" % (pprint.pformat(fnargs)))

def _combine_cwl_records(recs, fnargs):
    """Provide a list of nexted CWL records keyed by output key.
    """
    output_keys = _get_output_cwl_keys(fnargs)
    assert len(output_keys) == 1, output_keys
    return {output_keys[0]: recs}

def _collapse_to_cwl_record(samples, work_dir):
    """Convert nested samples from batches into a CWL record, based on input keys.
    """
    input_keys = sorted(list(set().union(*[d["cwl_keys"] for d in samples])), key=lambda x: (-len(x), tuple(x)))
    out = {}
    for key in input_keys:
        key_parts = key.split("__")
        vals = []
        cur = []
        for d in samples:
            vals.append(_to_cwl(tz.get_in(key_parts, d)))
            # Remove nested keys to avoid specifying multiple times
            cur.append(_dissoc_in(d, key_parts) if len(key_parts) > 1 else d)
        samples = cur
        out[key] = vals
    return out

def _to_cwl(val):
    """Convert a value into CWL formatted JSON, handling files and complex things.
    """
    # files where we list the entire directory as secondary files
    dir_targets = ["mainIndex"]
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
            cur_dir, cur_file = os.path.split(val["path"])
            if cur_file in dir_targets:
                for fname in os.listdir(cur_dir):
                    if fname != cur_file:
                        secondary.append({"class": "File", "path": os.path.join(cur_dir, fname)})
            if secondary:
                val["secondaryFiles"] = secondary
    elif isinstance(val, (list, tuple)):
        val = [_to_cwl(x) for x in val]
    elif isinstance(val, dict):
        # File representation with secondary files
        if "base" in val and "secondary" in val:
            out = {"class": "File", "path": val["base"]}
            secondary = [{"class": "File", "path": x} for x in val["secondary"]]
            if secondary:
                out["secondaryFiles"] = secondary
            val = out
        else:
            val = json.dumps(val, sort_keys=True, separators=(',', ':'))
    return val

def _dissoc_in(d, key_parts):
    if len(key_parts) > 1:
        d[key_parts[0]] = _dissoc_in(d[key_parts[0]], key_parts[1:])
    else:
        d.pop(key_parts[0], None)
    return d

def _update_nested(key, val, data):
    """Update the data object, avoiding over-writing with nested dictionaries.
    """
    if isinstance(val, dict):
        for sub_key, sub_val in val.items():
            data = _update_nested(key + [sub_key], sub_val, data)
    else:
        already_there = tz.get_in(key, data) is not None
        if already_there and val:
            raise ValueError("Duplicated key %s: %s and %s" % (key, val, tz.get_in(key, data)))
        if val or not already_there:
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
