"""Run distributed functions provided a name and json/YAML file with arguments.

Enables command line access and alternative interfaces to run specific
functionality within bcbio-nextgen.
"""
import collections
import contextlib
import csv
import json
import operator
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
        fnargs, parallel, out_keys = _world_from_cwl(args.name, fnargs[1:], work_dir)
        # Can remove this awkward Docker merge when we do not need custom GATK3 installs
        fnargs = config_utils.merge_resources(fnargs)
        argfile = os.path.join(work_dir, "cwl.output.json")
    else:
        parallel, out_keys = None, {}
    with utils.chdir(work_dir):
        with contextlib.closing(log.setup_local_logging(parallel={"wrapper": "runfn"})):
            try:
                out = fn(fnargs)
            except:
                logger.exception()
                raise
    if argfile:
        try:
            _write_out_argfile(argfile, out, fnargs, parallel, out_keys, work_dir)
        except:
            logger.exception()
            raise
        if argfile.endswith(".json"):
            _write_wdl_outputs(argfile, out_keys)

def _write_wdl_outputs(argfile, out_keys):
    """Write variables as WDL compatible output files.

    Writes individual files prefixed with 'wdl.output' that can be read
    by WDL standard library functions:

    https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#outputs
    """
    out_basename = "wdl.output.%s.txt"
    with open(argfile) as in_handle:
        outputs = json.load(in_handle)
    record_name, record_attrs = _get_record_attrs(out_keys)
    if record_name:
        recs = outputs[record_name]
        with open(out_basename % record_name, "w") as out_handle:
            writer = csv.writer(out_handle)
            if not isinstance(recs, (list, tuple)):
                recs = [recs]
            recs = list(utils.flatten(recs))
            keys = sorted(list(set(reduce(operator.add, [r.keys() for r in recs]))))
            writer.writerow(keys)
            for rec in recs:
                writer.writerow([_cwlvar_to_wdl(rec.get(k)) for k in keys])
    else:
        for key in out_keys:
            with open(out_basename % key, "w") as out_handle:
                vals = _cwlvar_to_wdl(outputs.get(key))
                if not isinstance(vals, (list, tuple)):
                    vals = [vals]
                for val in vals:
                    if isinstance(val, (list, tuple)):
                        val = "\t".join([str(x) for x in val])
                    out_handle.write(str(val) + "\n")

def _cwlvar_to_wdl(var):
    """Convert a CWL output object into a WDL output.

    This flattens files and other special CWL outputs that are
    plain strings in WDL.
    """
    if isinstance(var, (list, tuple)):
        return [_cwlvar_to_wdl(x) for x in var]
    elif isinstance(var, dict):
        assert var.get("class") == "File", var
        # XXX handle secondary files
        return var["path"]
    else:
        return var

def _write_out_argfile(argfile, out, fnargs, parallel, out_keys, work_dir):
    """Write output argfile, preparing a CWL ready JSON or YAML representation of the world.
    """
    with open(argfile, "w") as out_handle:
        if argfile.endswith(".json"):
            record_name, record_attrs = _get_record_attrs(out_keys)
            if record_name:
                if parallel in ["multi-batch"]:
                    recs = _nested_cwl_record(out, record_attrs)
                elif parallel in ["single-split", "multi-combined", "multi-parallel", "batch-single"]:
                    recs = [_collapse_to_cwl_record_single(utils.to_single_data(xs), record_attrs) for xs in out]
                else:
                    samples = [utils.to_single_data(xs) for xs in out]
                    recs = [_collapse_to_cwl_record(samples, record_attrs)]
                json.dump(_combine_cwl_records(recs, record_name, parallel),
                            out_handle, sort_keys=True, indent=4, separators=(', ', ': '))
            elif parallel in ["single-split", "multi-combined", "batch-split"]:
                json.dump(_convert_to_cwl_json([utils.to_single_data(xs) for xs in out], fnargs),
                            out_handle, sort_keys=True, indent=4, separators=(', ', ': '))
            else:
                json.dump(_convert_to_cwl_json(utils.to_single_data(utils.to_single_data(out)), fnargs),
                            out_handle, sort_keys=True, indent=4, separators=(', ', ': '))
        else:
            yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)

def _get_record_attrs(out_keys):
    """Check for records, a single key plus output attributes.
    """
    if len(out_keys) == 1:
        attr = out_keys.keys()[0]
        if out_keys[attr]:
            return attr, out_keys[attr]
    return None, None

def _add_resources(data, runtime):
    if "config" in data:
        data["config"]["resources"] = {"default": {"cores": int(runtime["cores"]),
                                                   "memory": "%sM" % int(float(runtime["ram"]) /
                                                                         float(runtime["cores"]))}}
        data["config"]["algorithm"]["num_cores"] = int(runtime["cores"])
    return data

def _world_from_cwl(fn_name, fnargs, work_dir):
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
        if key == "sentinel_parallel":
            parallel = val
            continue
        if key == "sentinel_runtime":
            runtime = dict(tz.partition(2, val.split(",")))
            continue
        if key == "sentinel_outputs":
            output_cwl_keys = _parse_output_keys(val)
            continue
        if key == "sentinel_inputs":
            input_order = collections.OrderedDict([x.split(":") for x in val.split(",")])
            continue
        else:
            assert key not in passed_keys, "Multiple keys should be handled via JSON records"
            passed_keys.append(key)
            key = key.split("__")
            data = _update_nested(key, _convert_value(val), data)
    if data:
        out.append(_finalize_cwl_in(data, work_dir, passed_keys, output_cwl_keys, runtime))

    # Read inputs from standard files instead of command line
    assert os.path.exists(os.path.join(work_dir, "cwl.inputs.json"))
    out = _read_from_cwlinput(os.path.join(work_dir, "cwl.inputs.json"), work_dir, runtime, parallel,
                              input_order, output_cwl_keys)

    if parallel in ["single-parallel", "single-merge", "multi-parallel", "multi-combined", "multi-batch",
                    "batch-split", "batch-parallel", "batch-merge", "batch-single"]:
        out = [out]
    else:
        assert len(out) == 1, "%s\n%s" % (pprint.pformat(out), pprint.pformat(fnargs))
    return out, parallel, output_cwl_keys

def _parse_output_keys(val):
    """Parse expected output keys from string, handling records.
    """
    out = {}
    for k in val.split(","):
        # record output
        if ":" in k:
            name, attrs = k.split(":")
            out[name] = attrs.split(";")
        else:
            out[k] = None
    return out

def _read_from_cwlinput(in_file, work_dir, runtime, parallel, input_order, output_cwl_keys):
    """Read data records from a JSON dump of inputs. Avoids command line flattening of records.
    """
    with open(in_file) as in_handle:
        inputs = json.load(in_handle)
    items_by_key = {}
    passed_keys = set([])
    for key, input_val in ((k, v) for (k, v) in inputs.items() if not k.startswith(("sentinel", "ignore"))):
        if key.endswith("_toolinput"):
            key = key.replace("_toolinput", "")
        if input_order[key] == "record":
            cur_keys, items = _read_cwl_record(input_val)
            passed_keys |= cur_keys
            items_by_key[key] = items
        else:
            items_by_key[tuple(key.split("__"))] = _cwlvar_to_wdl(input_val)
    prepped = _merge_cwlinputs(items_by_key, input_order, parallel)
    out = []
    for data in prepped:
        if isinstance(data, (list, tuple)):
            out.append([_finalize_cwl_in(utils.to_single_data(x), work_dir, list(passed_keys),
                                         output_cwl_keys, runtime) for x in data])
        else:
            out.append(_finalize_cwl_in(data, work_dir, list(passed_keys), output_cwl_keys, runtime))
    return out

def _merge_cwlinputs(items_by_key, input_order, parallel):
    """Merge multiple cwl records and inputs, handling multiple data items.

    Special cases:
    - Single record but multiple variables (merging arrayed jobs). Assign lists
      of variables to the record.
    """
    def item_count(x):
        return len(x) if isinstance(x, (list, tuple)) else 1
    var_items = set([item_count(items_by_key[tuple(k.split("__"))]) for (k, t) in input_order.items() if t == "var"])
    rec_items = set([item_count(items_by_key[k]) for (k, t) in input_order.items() if t == "record"])
    if var_items:
        num_items = var_items
        assert len(num_items) == 1, "Non-consistent variant data counts in CWL input:\n%s" % \
            (pprint.pformat(items_by_key))
        items_by_key, num_items = _nest_vars_in_rec(var_items, rec_items, input_order, items_by_key, parallel)
    else:
        num_items = rec_items
        assert len(num_items) == 1, "Non-consistent record data counts in CWL input:\n%s" % \
            (pprint.pformat(items_by_key))
    out = [{} for _ in range(num_items.pop())]
    for (cwl_key, cwl_type) in input_order.items():
        if cwl_type == "var":
            cwl_key = tuple(cwl_key.split("__"))
        cur_vals = items_by_key[cwl_key]
        for i, cur_val in enumerate(cur_vals):
            if isinstance(cwl_key, (list, tuple)):
                # nested batches with records
                if parallel.startswith("batch") and isinstance(out[i], (list, tuple)):
                    for j in range(len(out[i])):
                        out[i][j] = _update_nested(list(cwl_key), cur_val, out[i][j])
                else:
                    out[i] = _update_nested(list(cwl_key), cur_val, out[i])
            elif out[i] == {}:
                out[i] = cur_val
            else:
                assert isinstance(cur_val, dict), (cwl_key, cur_val)
                for k, v in cur_val.items():
                    out[i] = _update_nested([k], v, out[i], allow_overwriting=True)
    return out

def _nest_vars_in_rec(var_items, rec_items, input_order, items_by_key, parallel):
    """Nest multiple variable inputs into a single record or list of batch records.

    Custom CWL implementations extract and merge these.
    """
    num_items = var_items
    var_items = list(var_items)[0]
    if rec_items:
        rec_items = list(rec_items)[0]
        if ((rec_items == 1 and var_items > 1) or parallel.startswith("batch")):
            num_items = set([rec_items])
            for var_key in (k for (k, t) in input_order.items() if t != "record"):
                var_key = tuple(var_key.split("__"))
                items_by_key[var_key] = [items_by_key[var_key]] * rec_items
        else:
            assert var_items == rec_items, (var_items, rec_items)
    return items_by_key, num_items

def _expand_rec_to_vars(var_items, rec_items, input_order, items_by_key, parallel):
    """Expand record to apply to number of variants.

    Alternative approach to _nest_vars_in_rec
    to combining a single record with multiple variants.
    """
    num_items = var_items
    var_items = list(var_items)[0]
    if rec_items:
        for rec_key in (k for (k, t) in input_order.items() if t == "record"):
            rec_vals = items_by_key[rec_key]
            if len(rec_vals) == 1 and var_items > 1:
                items_by_key[rec_key] = rec_vals * var_items
            else:
                assert var_items == len(rec_vals), (var_items, rec_vals)
    return items_by_key, num_items

def _read_cwl_record(rec):
    """Read CWL records, handling multiple nesting and batching cases.
    """
    keys = set([])
    out = []
    if isinstance(rec, dict):
        is_batched = all([isinstance(v, (list, tuple)) for v in rec.values()])
        cur = [{} for _ in range(len(rec.values()[0]) if is_batched else 1)]
        for k in rec.keys():
            keys.add(k)
            val = rec[k]
            val = val if is_batched else [val]
            for i, v in enumerate(val):
                v = _cwlvar_to_wdl(v)
                cur[i] = _update_nested(k.split("__"), v, cur[i])
        if is_batched:
            out.append(cur)
        else:
            assert len(cur) == 1
            out.append(cur[0])
    else:
        assert isinstance(rec, (list, tuple))
        for sub_rec in rec:
            sub_keys, sub_out = _read_cwl_record(sub_rec)
            keys |= sub_keys
            out.append(sub_out)
    return keys, out

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
    if isinstance(val, (list, tuple)):
        return [_convert_value(x) for x in val]
    elif val is None:
        return val
    elif _is_number(val, int):
        return int(val)
    elif _is_number(val, float):
        return float(val)
    elif val.find(";;") >= 0:
        return [_convert_value(v) for v in val.split(";;")]
    elif val.startswith(("{", "[")):
        # Can get ugly JSON output from CWL with unicode and ' instead of "
        # This tries to fix it so parsed correctly by json loader
        return json.loads(val.replace("u'", "'").replace("'", '"'))
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
        keys = []
        for key in outvar.split("__"):
            try:
                key = int(key)
            except ValueError:
                pass
            keys.append(key)
        if isinstance(data, dict):
            out[outvar] = _to_cwl(tz.get_in(keys, data))
        else:
            out[outvar] = [_to_cwl(tz.get_in(keys, x)) for x in data]
    return out

def _get_output_cwl_keys(fnargs):
    """Retrieve output_cwl_keys from potentially nested input arguments.
    """
    for d in utils.flatten(fnargs):
        if isinstance(d, dict) and d.get("output_cwl_keys"):
            return d["output_cwl_keys"]
    raise ValueError("Did not find output_cwl_keys in %s" % (pprint.pformat(fnargs)))

def _combine_cwl_records(recs, record_name, parallel):
    """Provide a list of nexted CWL records keyed by output key.

    Handles batches, where we return a list of records, and single items
    where we return one record.
    """
    if parallel not in ["multi-batch", "single-split", "multi-combined", "batch-single"]:
        assert len(recs) == 1, pprint.pformat(recs)
        return {record_name: recs[0]}
    else:
        return {record_name: recs}

def _collapse_to_cwl_record_single(data, want_attrs):
    """Convert a single sample into a CWL record.
    """
    out = {}
    for key in want_attrs:
        key_parts = key.split("__")
        out[key] = _to_cwl(tz.get_in(key_parts, data))
    return out

def _nested_cwl_record(xs, want_attrs):
    """Convert arbitrarily nested samples into a nested list of dictionaries.

    nests only at the record level, rather than within records. For batching
    a top level list is all of the batches and sub-lists are samples within the
    batch.
    """
    if isinstance(xs, (list, tuple)):
        return [_nested_cwl_record(x, want_attrs) for x in xs]
    else:
        assert isinstance(xs, dict), pprint.pformat(xs)
        return _collapse_to_cwl_record_single(xs, want_attrs)

def _collapse_to_cwl_record(samples, want_attrs):
    """Convert nested samples from batches into a CWL record, based on input keys.
    """
    input_keys = sorted(list(set().union(*[d["cwl_keys"] for d in samples])), key=lambda x: (-len(x), tuple(x)))
    out = {}
    for key in input_keys:
        if key in want_attrs:
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
    # aligner and database indices where we list the entire directory as secondary files
    dir_targets = ("mainIndex", ".alt", ".amb", ".ann", ".bwt", ".pac", ".sa", ".ebwt", ".bt2",
                   "Genome", "GenomeIndex", "GenomeIndexHash", "OverflowTable")
    if isinstance(val, basestring):
        if os.path.exists(val) and os.path.isfile(val):
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
            # Handle relative paths
            if not cur_dir:
                cur_dir = os.getcwd()
            if cur_file.endswith(dir_targets):
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

def _update_nested(key, val, data, allow_overwriting=False):
    """Update the data object, avoiding over-writing with nested dictionaries.
    """
    if isinstance(val, dict):
        for sub_key, sub_val in val.items():
            data = _update_nested(key + [sub_key], sub_val, data)
    else:
        already_there = tz.get_in(key, data) is not None
        if already_there and val:
            if not allow_overwriting:
                raise ValueError("Duplicated key %s: %s and %s" % (key, val, tz.get_in(key, data)))
            else:
                already_there = False
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
