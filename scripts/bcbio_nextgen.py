#!/usr/bin/env python -Es
"""Run an automated analysis pipeline for high throughput sequencing data.

Handles runs in local or distributed mode based on the command line or
configured parameters.

The <config file> is a global YAML configuration file specifying details
about the system. An example configuration file is in 'config/bcbio_system.yaml'.
This is optional for automated installations.

<fc_dir> is an optional parameter specifying a directory of Illumina output
or fastq files to process. If configured to connect to a Galaxy LIMS system,
this can retrieve run information directly from Galaxy for processing.

<YAML run information> is an optional file that specifies details about the
flowcell lanes, instead of retrieving it from Galaxy. An example
configuration file is located in 'config/bcbio_sample.yaml' This allows running
on files in arbitrary locations with no connection to Galaxy required.

Usage:
  bcbio_nextgen.py <config_file> [<fc_dir>] [<run_info_yaml>]
     -t type of parallelization to use:
          - local: Non-distributed, possibly multiple if n > 1 (default)
          - ipython: IPython distributed processing
     -n total number of processes to use
     -s scheduler for ipython parallelization (lsf, sge, slurm, torque, pbspro)
     -q queue to submit jobs for ipython parallelization
"""
import os
import argparse
import sys

from bcbio import install, utils, workflow
from bcbio.illumina import machine
from bcbio.distributed import runfn, clargs
from bcbio.pipeline.main import run_main
from bcbio.server import main as server_main
from bcbio.graph import graph
from bcbio.provenance import programs
from bcbio.pipeline import version

def main(**kwargs):
    run_main(**kwargs)

def parse_cl_args(in_args):
    """Parse input commandline arguments, handling multiple cases.

    Returns the main config file and set of kwargs.
    """
    sub_cmds = {"upgrade": install.add_subparser,
                "server": server_main.add_subparser,
                "runfn": runfn.add_subparser,
                "graph": graph.add_subparser,
                "version": programs.add_subparser,
                "sequencer": machine.add_subparser}
    description = "Community developed high throughput sequencing analysis."
    parser = argparse.ArgumentParser(description=description)
    sub_cmd = None
    if len(in_args) > 0 and in_args[0] in sub_cmds:
        subparser_help = "bcbio-nextgen supplemental commands"
        subparsers = parser.add_subparsers(help=subparser_help)
        sub_cmds[in_args[0]](subparsers)
        sub_cmd = in_args[0]
    else:
        parser.add_argument("global_config", nargs="?",
                            help=("Global YAML configuration file specifying "
                                  "details about the system (optional, "
                                  "defaults to installed bcbio_system.yaml)"))
        parser.add_argument("fc_dir", nargs="?",
                            help=("A directory of Illumina output or fastq "
                                  "files to process (optional)"))
        parser.add_argument("run_config", nargs="*",
                            help=("YAML file with details about samples to "
                                  "process (required, unless using Galaxy "
                                  "LIMS as input)")),
        parser.add_argument("-n", "--numcores", type=int, default=1,
                            help="Total cores to use for processing")
        parser.add_argument("-t", "--paralleltype",
                            choices=["local", "ipython"],
                            default="local", help="Approach to parallelization")
        parser.add_argument("-s", "--scheduler",
                            choices=["lsf", "sge", "torque", "slurm", "pbspro"],
                            help="Scheduler to use for ipython parallel")
        parser.add_argument("-q", "--queue",
                            help=("Scheduler queue to run jobs on, for "
                                  "ipython parallel"))
        parser.add_argument("-r", "--resources",
                            help=("Cluster specific resources specifications. "
                                  "Can be specified multiple times.\n"
                                  "Supports SGE, Torque, LSF and SLURM "
                                  "parameters."), default=[], action="append")
        parser.add_argument("--timeout", default=15, type=int,
                            help=("Number of minutes before cluster startup "
                                  "times out. Defaults to 15"))
        parser.add_argument("--retries", default=0, type=int,
                            help=("Number of retries of failed tasks during "
                                  "distributed processing. Default 0 "
                                  "(no retries)"))
        parser.add_argument("-p", "--tag",
                            help="Tag name to label jobs on the cluster",
                            default="")
        parser.add_argument("-w", "--workflow",
                            help=("Run a workflow with the given commandline "
                                  "arguments"))
        parser.add_argument("--workdir", default=os.getcwd(),
                            help=("Directory to process in. Defaults to "
                                  "current working directory"))
        parser.add_argument("-v", "--version", help="Print current version",
                            action="store_true")
        # Hidden arguments passed downstream
        parser.add_argument("--only-metadata", help=argparse.SUPPRESS, action="store_true", default=False)
    args = parser.parse_args(in_args)
    if hasattr(args, "workdir"):
        args.workdir = utils.safe_makedir(os.path.abspath(args.workdir))
    if hasattr(args, "global_config"):
        error_msg = _sanity_check_args(args)
        if error_msg:
            parser.error(error_msg)
        kwargs = {"parallel": clargs.to_parallel(args),
                  "workflow": args.workflow,
                  "workdir": args.workdir}
        kwargs = _add_inputs_to_kwargs(args, kwargs, parser)
        error_msg = _sanity_check_kwargs(kwargs)
        if error_msg:
            parser.error(error_msg)
    else:
        assert sub_cmd is not None
        kwargs = {"args": args,
                  "config_file": None,
                  sub_cmd: True}
    return kwargs

def _sanity_check_args(args):
    """Ensure dependent arguments are correctly specified
    """
    if "scheduler" in args and "queue" in args:
        if args.scheduler and not args.queue:
            if args.scheduler != "sge":
                return "IPython parallel scheduler (-s) specified. This also requires a queue (-q)."
        elif args.queue and not args.scheduler:
            return "IPython parallel queue (-q) supplied. This also requires a scheduler (-s)."
        elif args.paralleltype == "ipython" and (not args.queue or not args.scheduler):
            return "IPython parallel requires queue (-q) and scheduler (-s) arguments."

def _sanity_check_kwargs(args):
    """Sanity check after setting up input arguments, handling back compatibility
    """
    if not args.get("workflow") and not args.get("run_info_yaml"):
        return ("Require a sample YAML file describing inputs: "
                "https://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html")

def _add_inputs_to_kwargs(args, kwargs, parser):
    """Convert input system config, flow cell directory and sample yaml to kwargs.

    Handles back compatibility with previous commandlines while allowing flexible
    specification of input parameters.
    """
    inputs = [x for x in [args.global_config, args.fc_dir] + args.run_config
              if x is not None]
    global_config = "bcbio_system.yaml"  # default configuration if not specified
    if kwargs.get("workflow", "") == "template":
        if args.only_metadata:
            inputs.append("--only-metadata")
        kwargs["inputs"] = inputs
        return kwargs
    elif len(inputs) == 1:
        if os.path.isfile(inputs[0]):
            fc_dir = None
            run_info_yaml = inputs[0]
        else:
            fc_dir = inputs[0]
            run_info_yaml = None
    elif len(inputs) == 2:
        if os.path.isfile(inputs[0]):
            global_config = inputs[0]
            if os.path.isfile(inputs[1]):
                fc_dir = None
                run_info_yaml = inputs[1]
            else:
                fc_dir = inputs[1]
                run_info_yaml = None
        else:
            fc_dir, run_info_yaml = inputs
    elif len(inputs) == 3:
        global_config, fc_dir, run_info_yaml = inputs
    elif args.version:
        print version.__version__
        sys.exit()
    else:
        print "Incorrect input arguments", inputs
        parser.print_help()
        sys.exit()
    if fc_dir:
        fc_dir = os.path.abspath(fc_dir)
    if run_info_yaml:
        run_info_yaml = os.path.abspath(run_info_yaml)
    if kwargs.get("workflow"):
        kwargs["inputs"] = inputs
    kwargs["config_file"] = global_config
    kwargs["fc_dir"] = fc_dir
    kwargs["run_info_yaml"] = run_info_yaml
    return kwargs

if __name__ == "__main__":
    kwargs = parse_cl_args(sys.argv[1:])
    if "upgrade" in kwargs and kwargs["upgrade"]:
        install.upgrade_bcbio(kwargs["args"])
    elif "server" in kwargs and kwargs["server"]:
        server_main.start(kwargs["args"])
    elif "runfn" in kwargs and kwargs["runfn"]:
        runfn.process(kwargs["args"])
    elif "graph" in kwargs and kwargs["graph"]:
        graph.bootstrap(kwargs["args"])
    elif "version" in kwargs and kwargs["version"]:
        programs.write_versions({"work": kwargs["args"].workdir})
    elif "sequencer" in kwargs and kwargs["sequencer"]:
        machine.check_and_postprocess(kwargs["args"])
    else:
        if kwargs.get("workflow"):
            setup_info = workflow.setup(kwargs["workflow"], kwargs.pop("inputs"))
            if setup_info is None:  # no automated run after setup
                sys.exit(0)
            workdir, new_kwargs = setup_info
            os.chdir(workdir)
            kwargs.update(new_kwargs)
        main(**kwargs)

