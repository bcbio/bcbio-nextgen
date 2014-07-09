"""Run distributed functions provided a name and json/YAML file with arguments.

Enables command line access and alternative interfaces to run specific
functionality within bcbio-nextgen.
"""
import os

import yaml

from bcbio import log, utils
from bcbio.distributed import multitasks
from bcbio.pipeline import config_utils

def process(args):
    """Run the function in args.name given arguments in args.argfile.
    """
    try:
        fn = getattr(multitasks, args.name)
    except AttributeError:
        raise AttributeError("Did not find exposed function in bcbio.distributed.multitasks named '%s'" % args.name)
    with open(args.argfile) as in_handle:
        fnargs = yaml.safe_load(in_handle)
    fnargs = config_utils.merge_resources(fnargs)
    work_dir = os.path.dirname(args.argfile)
    with utils.chdir(work_dir):
        log.setup_local_logging(parallel={"wrapper": "runfn"})
        out = fn(fnargs)
    out_file = "%s-out%s" % os.path.splitext(args.argfile)
    with open(out_file, "w") as out_handle:
        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)

def add_subparser(subparsers):
    parser = subparsers.add_parser("runfn", help=("Run a specific bcbio-nextgen function."
                                                  "Intended for distributed use."))
    parser.add_argument("name", help="Name of the function to run")
    parser.add_argument("argfile", help="JSON file with arguments to the function")
