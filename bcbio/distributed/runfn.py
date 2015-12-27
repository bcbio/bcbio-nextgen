"""Run distributed functions provided a name and json/YAML file with arguments.

Enables command line access and alternative interfaces to run specific
functionality within bcbio-nextgen.
"""
import json
import os

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
        fnargs = _world_from_cwl(fnargs[1:], work_dir)
        argfile = os.path.join(work_dir, "cwl-%s-world.json" % args.name)
    with utils.chdir(work_dir):
        log.setup_local_logging(parallel={"wrapper": "runfn"})
        out = fn(fnargs)
    if argfile:
        with open(argfile, "w") as out_handle:
            if argfile.endswith(".json"):
                json.dump(_remove_work_dir(out[0][0], work_dir + "/"), out_handle)
            else:
                yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)

def _world_from_cwl(fnargs, work_dir):
    """Reconstitute a bcbio world data object from flattened CWL-compatible inputs.

    Converts the flat CWL representation into a nested bcbio world dictionary.
    """
    data = {}
    for fnarg in fnargs:
        key, val = fnarg.split("=")
        key = key.split("__")
        if val.startswith(("{", "[")):
            val = json.loads(val)
        elif val.find(";;") >= 0:
            val = val.split(";;")
        data = tz.update_in(data, key, lambda x: val)
    data["dirs"] = {"work": work_dir}
    # XXX Determine cores and other resources from CWL
    data["config"]["resources"] = {}
    data = run_info.normalize_world(data)
    return [data]

def _remove_work_dir(orig, work_dir):
    """Remove work directory from any file paths to make them relative.

    CWL needs to work off relative paths since files change on reloads.
    """
    def startswith_work_dir(x):
        return x and isinstance(x, basestring) and x.startswith(work_dir)
    out = {}
    for key, val in orig.items():
        if isinstance(val, dict):
            out[key] = _remove_work_dir(val, work_dir)
        elif isinstance(val, (list, tuple)):
            out[key] = [x.replace(work_dir, "") if startswith_work_dir(x) else x for x in val]
        else:
            out[key] = val.replace(work_dir, "") if startswith_work_dir(val) else val
    return out

def add_subparser(subparsers):
    parser = subparsers.add_parser("runfn", help=("Run a specific bcbio-nextgen function."
                                                  "Intended for distributed use."))
    parser.add_argument("name", help="Name of the function to run")
    parser.add_argument("argfile", help="JSON file with arguments to the function")
    parser.add_argument("moreargs", nargs="*", help="Additional arguments to pass in the case of raw input")
    parser.add_argument("--raw", action="store_true", default=False,
                        help="Treat the inputs as raw file arguments to the function, instead of parsing them.")
    parser.add_argument("-o", "--outfile", help="Output file to write, defaults to inputfile-out.yaml")
