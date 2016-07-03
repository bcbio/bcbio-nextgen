"""Run bcbio generated CWL with a supported tool.

Handles wrapping and integrating with multiple tools making it easier
to run bcbio in a standard way in many environments.
"""
import glob
import os
import subprocess
import sys

from bcbio.provenance import do

def _get_main_and_json(directory):
    main_cwl = glob.glob(os.path.join(directory, "main-*.cwl"))
    assert len(main_cwl) == 1, "Did not find main CWL in %s" % directory
    main_json = glob.glob(os.path.join(directory, "main-*-samples.json"))
    assert len(main_json) == 1, "Did not find main json in %s" % directory
    return main_cwl[0], main_json[0]

def _run_cwltool(directory, args):
    """Run with cwltool -- reference implementation.
    """
    main_file, json_file = _get_main_and_json(directory)
    cmd = [os.path.join(os.path.dirname(sys.executable), "cwltool")] + args + \
          [main_file, json_file]
    subprocess.check_call(cmd)
    #do.run(cmd, "Run bcbio CWL with cwltool")

def _run_arvados(directory, args):
    """Run CWL on Aravdos.
    """
    pass

def _run_toil(directory, args):
    """Run CWL with Toil.
    """
    pass

_TOOLS = {"cwltool": _run_cwltool,
          "arvados": _run_arvados,
          "toil": _run_toil}

def run(args):
    _TOOLS[args.tool](args.directory, args.toolargs)
