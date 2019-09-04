"""Run distributed functions provided a name and json/YAML file with arguments.

Enables command line access and alternative interfaces to run specific
functionality within bcbio-nextgen.
"""
import collections
import contextlib
import json
import os
import pprint
import shutil

import six
import toolz as tz
import yaml

from bcbio import log, utils, setpath, utils
from bcbio.log import logger
from bcbio.cwl import cwlutils
from bcbio.distributed import multitasks
from bcbio.pipeline import config_utils, run_info

def process(args):
    """Run the function in args.name given arguments in args.argfile.
    """
    # Set environment to standard to use periods for decimals and avoid localization
    locale_to_use = utils.get_locale()
    os.environ["LC_ALL"] = locale_to_use
    os.environ["LC"] = locale_to_use
    os.environ["LANG"] = locale_to_use
    setpath.prepend_bcbiopath()
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
        fnargs, parallel, out_keys, input_files = _world_from_cwl(args.name, fnargs[1:], work_dir)
        # Can remove this awkward Docker merge when we do not need custom GATK3 installs
        fnargs = config_utils.merge_resources(fnargs)
        argfile = os.path.join(work_dir, "cwl.output.json")
    else:
        parallel, out_keys, input_files = None, {}, []
    with utils.chdir(work_dir):
        with contextlib.closing(log.setup_local_logging(parallel={"wrapper": "runfn"})):
            try:
                out = fn(*fnargs)
            except:
                logger.exception()
                raise
            finally:
                # Clean up any copied and unpacked workflow inputs, avoiding extra disk usage
                wf_input_dir = os.path.join(work_dir, "wf-inputs")
                if os.path.exists(wf_input_dir) and os.path.isdir(wf_input_dir):
                    shutil.rmtree(wf_input_dir)
    if argfile:
        try:
            _write_out_argfile(argfile, out, fnargs, parallel, out_keys, input_files, work_dir)
        except:
            logger.exception()
            raise

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
        return var.get("path") or var["value"]
    else:
        return var

def _write_out_argfile(argfile, out, fnargs, parallel, out_keys, input_files, work_dir):
    """Write output argfile, preparing a CWL ready JSON or YAML representation of the world.
    """
    with open(argfile, "w") as out_handle:
        if argfile.endswith(".json"):
            record_name, record_attrs = _get_record_attrs(out_keys)
            if record_name:
                if parallel in ["multi-batch"]:
                    recs = _nested_cwl_record(out, record_attrs, input_files)
                elif parallel in ["single-split", "multi-combined", "multi-parallel", "batch-single",
                                  "single-single"]:
                    recs = [_collapse_to_cwl_record_single(utils.to_single_data(xs), record_attrs, input_files)
                            for xs in out]
                else:
                    samples = [utils.to_single_data(xs) for xs in out]
                    recs = [_collapse_to_cwl_record(samples, record_attrs, input_files)]
                json.dump(_combine_cwl_records(recs, record_name, parallel),
                            out_handle, sort_keys=True, indent=4, separators=(', ', ': '))
            elif parallel in ["single-split", "multi-combined", "batch-split"]:
                json.dump(_convert_to_cwl_json([utils.to_single_data(xs) for xs in out], fnargs, input_files),
                            out_handle, sort_keys=True, indent=4, separators=(', ', ': '))
            else:
                json.dump(_convert_to_cwl_json(utils.to_single_data(utils.to_single_data(out)), fnargs, input_files),
                            out_handle, sort_keys=True, indent=4, separators=(', ', ': '))
        else:
            yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)

def _get_record_attrs(out_keys):
    """Check for records, a single key plus output attributes.
    """
    if len(out_keys) == 1:
        attr = list(out_keys.keys())[0]
        if out_keys[attr]:
            return attr, out_keys[attr]
    return None, None

def _add_resources(data, runtime):
    """Merge input resources with current CWL runtime parameters.
    """
    if "config" not in data:
        data["config"] = {}
    # Convert input resources, which may be a JSON string
    resources = data.get("resources", {}) or {}
    if isinstance(resources, six.string_types) and resources.startswith(("{", "[")):
        resources = json.loads(resources)
        data["resources"] = resources
    assert isinstance(resources, dict), (resources, data)
    data["config"]["resources"] = resources
    # Add in memory and core usage from CWL
    memory = int(float(runtime["ram"]) / float(runtime["cores"]))
    data["config"]["resources"].update({"default": {"cores": int(runtime["cores"]),
                                                    "memory": "%sM" % memory,
                                                    "jvm_opts": ["-Xms%sm" % min(1000, memory // 2),
                                                                    "-Xmx%sm" % memory]}})
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
    out, input_files = _read_from_cwlinput(os.path.join(work_dir, "cwl.inputs.json"), work_dir, runtime, parallel,
                                           input_order, output_cwl_keys)

    if parallel in ["single-parallel", "single-merge", "multi-parallel", "multi-combined", "multi-batch",
                    "batch-split", "batch-parallel", "batch-merge", "batch-single"]:
        out = [out]
    else:
        assert len(out) == 1, "%s\n%s" % (pprint.pformat(out), pprint.pformat(fnargs))
    return out, parallel, output_cwl_keys, input_files

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

def _find_input_files(var, out):
    """Find input files within the given CWL object.
    """
    if isinstance(var, (list, tuple)):
        for x in var:
            out = _find_input_files(x, out)
    elif isinstance(var, dict):
        if var.get("class") == "File":
            out.append(var["path"])
            out = _find_input_files(var.get("secondaryFiles", []), out)
        for key, val in var.items():
            out = _find_input_files(val, out)
    return out

def _read_from_cwlinput(in_file, work_dir, runtime, parallel, input_order, output_cwl_keys):
    """Read data records from a JSON dump of inputs. Avoids command line flattening of records.
    """
    with open(in_file) as in_handle:
        inputs = json.load(in_handle)
    items_by_key = {}
    input_files = []
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
        input_files = _find_input_files(input_val, input_files)
    prepped = _merge_cwlinputs(items_by_key, input_order, parallel)
    out = []
    for data in prepped:
        if isinstance(data, (list, tuple)):
            out.append([_finalize_cwl_in(utils.to_single_data(x), work_dir, list(passed_keys),
                                         output_cwl_keys, runtime) for x in data])
        else:
            out.append(_finalize_cwl_in(data, work_dir, list(passed_keys), output_cwl_keys, runtime))
    return out, input_files

def _is_nested_item(x):
    return isinstance(x, (list, tuple))

def _maybe_nest_bare_single(items_by_key, parallel):
    """Nest single inputs to avoid confusing single items and lists like files.
    """
    if (parallel == "multi-parallel" and
          (sum([1 for x in items_by_key.values() if not _is_nested_item(x)]) >=
           sum([1 for x in items_by_key.values() if _is_nested_item(x)]))):
        out = {}
        for k, v in items_by_key.items():
            out[k] = [v]
        return out
    else:
        return items_by_key

def _item_count(x):
    return len(x) if _is_nested_item(x) else 1

def _is_nested_single(v, target):
    return target > 1 and _is_nested_item(v) and _item_count(v) == 1 and _item_count(v[0]) == target

def _check_for_single_nested(target, items_by_key, input_order):
    """Check for single nested inputs that match our target count and unnest.

    Handles complex var inputs where some have an extra layer of nesting.
    """
    out = utils.deepish_copy(items_by_key)
    for (k, t) in input_order.items():
        if t == "var":
            v = items_by_key[tuple(k.split("__"))]
            if _is_nested_single(v, target):
                out[tuple(k.split("__"))] = v[0]
    return out

def _concat_records(items_by_key, input_order):
    """Concatenate records into a single key to avoid merging.

    Handles heterogeneous records that will then be sorted out in
    the processing fuction.
    """
    all_records = []
    for (k, t) in input_order.items():
        if t == "record":
            all_records.append(k)
    out_items_by_key = utils.deepish_copy(items_by_key)
    out_input_order = utils.deepish_copy(input_order)
    if len(all_records) > 1:
        final_k = all_records[0]
        final_v = items_by_key[final_k]
        for k in all_records[1:]:
            final_v += items_by_key[k]
            del out_items_by_key[k]
            del out_input_order[k]
        out_items_by_key[final_k] = final_v
    return out_items_by_key, out_input_order

def _merge_cwlinputs(items_by_key, input_order, parallel):
    """Merge multiple cwl records and inputs, handling multiple data items.

    Special cases:
    - Single record but multiple variables (merging arrayed jobs). Assign lists
      of variables to the record.
    """
    items_by_key = _maybe_nest_bare_single(items_by_key, parallel)
    if parallel == "multi-combined":
        items_by_key, input_order = _concat_records(items_by_key, input_order)
    var_items = set([_item_count(items_by_key[tuple(k.split("__"))])
                     for (k, t) in input_order.items() if t == "var"])
    rec_items = set([_item_count(items_by_key[k]) for (k, t) in input_order.items() if t == "record"])
    if var_items:
        num_items = var_items
        if len(num_items) == 2 and 1 in num_items:
            num_items.remove(1)
            items_by_key_test = _check_for_single_nested(num_items.pop(), items_by_key, input_order)
            var_items = set([_item_count(items_by_key_test[tuple(k.split("__"))])
                             for (k, t) in input_order.items() if t == "var"])
            num_items = var_items
        assert len(num_items) == 1, "Non-consistent variable data counts in CWL input:\n%s" % \
            (pprint.pformat(items_by_key))
        items_by_key, num_items = _nest_vars_in_rec(var_items, rec_items, input_order, items_by_key, parallel)
    else:
        num_items = rec_items
        assert len(num_items) == 1, "Non-consistent record data counts in CWL input:\n%s" % \
            (pprint.pformat(items_by_key))
    target_items = num_items.pop()
    out = [{} for _ in range(target_items)]
    for (cwl_key, cwl_type) in input_order.items():
        if cwl_type == "var":
            cwl_key = tuple(cwl_key.split("__"))
        cur_vals = items_by_key[cwl_key]
        if _is_nested_single(cur_vals, target_items):
            cur_vals = [[x] for x in cur_vals[0]]
        for i, cur_val in enumerate(cur_vals):
            if isinstance(cwl_key, (list, tuple)):
                # nested batches with records
                if (parallel.startswith(("batch", "multi-parallel")) and
                      isinstance(out[i], (list, tuple))):
                    for j in range(len(out[i])):
                        out[i][j] = _update_nested(list(cwl_key), cur_val, out[i][j], allow_overwriting=True)
                else:
                    out[i] = _update_nested(list(cwl_key), cur_val, out[i], allow_overwriting=True)
            elif out[i] == {}:
                out[i] = cur_val
            else:
                # Handle single non-batched records
                if isinstance(cur_val, (list, tuple)) and len(cur_val) == 1:
                    cur_val = cur_val[0]
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
    data = cwlutils.normalize_missing(data)
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

def _convert_to_cwl_json(data, fnargs, input_files):
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
            out[outvar] = _to_cwl(tz.get_in(keys, data), input_files)
        else:
            out[outvar] = [_to_cwl(tz.get_in(keys, x), input_files) for x in data]
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

def _collapse_to_cwl_record_single(data, want_attrs, input_files):
    """Convert a single sample into a CWL record.
    """
    out = {}
    for key in want_attrs:
        key_parts = key.split("__")
        out[key] = _to_cwl(tz.get_in(key_parts, data), input_files)
    return out

def _nested_cwl_record(xs, want_attrs, input_files):
    """Convert arbitrarily nested samples into a nested list of dictionaries.

    nests only at the record level, rather than within records. For batching
    a top level list is all of the batches and sub-lists are samples within the
    batch.
    """
    if isinstance(xs, (list, tuple)):
        return [_nested_cwl_record(x, want_attrs, input_files) for x in xs]
    else:
        assert isinstance(xs, dict), pprint.pformat(xs)
        return _collapse_to_cwl_record_single(xs, want_attrs, input_files)

def _collapse_to_cwl_record(samples, want_attrs, input_files):
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
                vals.append(_to_cwl(tz.get_in(key_parts, d), input_files))
                # Remove nested keys to avoid specifying multiple times
                cur.append(_dissoc_in(d, key_parts) if len(key_parts) > 1 else d)
            samples = cur
            out[key] = vals
    return out

def _file_and_exists(val, input_files):
    """Check if an input is a file and exists.

    Checks both locally (staged) and from input files (re-passed but never localized).
    """
    return ((os.path.exists(val) and os.path.isfile(val)) or
            val in input_files)

def _to_cwl(val, input_files):
    """Convert a value into CWL formatted JSON, handling files and complex things.
    """
    if isinstance(val, six.string_types):
        if _file_and_exists(val, input_files):
            val = {"class": "File", "path": val}
            secondary = []
            for idx in [".bai", ".tbi", ".gbi", ".fai", ".crai", ".db"]:
                idx_file = val["path"] + idx
                if _file_and_exists(idx_file, input_files):
                    secondary.append({"class": "File", "path": idx_file})
            for idx in [".dict"]:
                idx_file = os.path.splitext(val["path"])[0] + idx
                if _file_and_exists(idx_file, input_files):
                    secondary.append({"class": "File", "path": idx_file})
            cur_dir, cur_file = os.path.split(val["path"])
            # Handle relative paths
            if not cur_dir:
                cur_dir = os.getcwd()
            if cur_file.endswith(cwlutils.DIR_TARGETS):
                if os.path.exists(cur_dir):
                    for fname in os.listdir(cur_dir):
                        if fname != cur_file and not os.path.isdir(os.path.join(cur_dir, fname))\
                                and fname != 'sbg.worker.log':
                            secondary.append({"class": "File", "path": os.path.join(cur_dir, fname)})
                else:
                    for f in input_files:
                        if f.startswith(cur_dir) and f != cur_file and not os.path.isdir(f):
                            secondary.append({"class": "File", "path": f})
            if secondary:
                val["secondaryFiles"] = _remove_duplicate_files(secondary)
    elif isinstance(val, (list, tuple)):
        val = [_to_cwl(x, input_files) for x in val]
    elif isinstance(val, dict):
        # File representation with secondary files
        if "base" in val and "secondary" in val:
            out = {"class": "File", "path": val["base"]}
            secondary = [{"class": "File", "path": x} for x in val["secondary"] if not os.path.isdir(x)]
            if secondary:
                out["secondaryFiles"] = _remove_duplicate_files(secondary)
            val = out
        else:
            val = json.dumps(val, sort_keys=True, separators=(',', ':'))
    return val

def _remove_duplicate_files(xs):
    """Remove files specified multiple times in a list.
    """
    seen = set([])
    out = []
    for x in xs:
        if x["path"] not in seen:
            out.append(x)
            seen.add(x["path"])
    return out

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
            data = _update_nested(key + [sub_key], sub_val, data, allow_overwriting=allow_overwriting)
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
