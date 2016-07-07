"""Run bcbio generated CWL with a supported tool.

Handles wrapping and integrating with multiple tools making it easier
to run bcbio in a standard way in many environments.
"""
import glob
import os
import subprocess

from bcbio import utils

def _get_main_and_json(directory):
    """Retrieve the main CWL and sample JSON files from a bcbio generated directory.
    """
    directory = os.path.normpath(os.path.abspath(directory))
    main_cwl = glob.glob(os.path.join(directory, "main-*.cwl"))
    assert len(main_cwl) == 1, "Did not find main CWL in %s" % directory
    main_json = glob.glob(os.path.join(directory, "main-*-samples.json"))
    assert len(main_json) == 1, "Did not find main json in %s" % directory
    return main_cwl[0], main_json[0]

def _run_tool(cmd):
    if isinstance(cmd, (list, tuple)):
        cmd = " ".join([str(x) for x in cmd])
    cmd = utils.local_path_export() + cmd
    subprocess.check_call(cmd, shell=True)

def _run_cwltool(args):
    """Run with cwltool -- reference implementation.
    """
    main_file, json_file = _get_main_and_json(args.directory)
    if args.no_container:
        flags = ["--no-container", "--preserve-environment", "PATH", "HOME"]
    else:
        flags = []
    cmd = ["cwltool"] + flags + args.toolargs + [main_file, json_file]
    _run_tool(cmd)

def _run_arvados(args):
    """Run CWL on Aravdos.
    """
    assert not args.no_container, "No container is not available with Arvados runs"
    assert "ARVADOS_API_TOKEN" in os.environ and "ARVADOS_API_HOST" in os.environ, \
        "Need to set ARVADOS_API_TOKEN and ARVADOS_API_HOST in environment to run"
    main_file, json_file = _get_main_and_json(args.directory)
    flags = ["--local", "--enable-reuse"]
    cmd = ["arvados-cwl-runner"] + flags + args.toolargs + [main_file, json_file]
    _run_tool(cmd)

def _run_toil(args):
    """Run CWL with Toil.
    """
    main_file, json_file = _get_main_and_json(args.directory)
    work_dir = utils.safe_makedir(os.path.join(os.getcwd(), "cwltoil_work"))
    log_file = os.path.join(work_dir, "cwltoil.log")
    jobstore = os.path.join(work_dir, "cwltoil_jobstore")
    flags = ["--jobStore", "file:%s" % jobstore, "--logFile", log_file, "--workDir", work_dir]
    if os.path.exists(jobstore):
        flags += ["--restart"]
    if args.no_container:
        flags += ["--no-container", "--preserve-environment", "PATH", "HOME"]
    cmd = ["cwltoil"] + flags + args.toolargs + [main_file, json_file]
    with utils.chdir(work_dir):
        _run_tool(cmd)

_TOOLS = {"cwltool": _run_cwltool,
          "arvados": _run_arvados,
          "toil": _run_toil}

def run(args):
    _TOOLS[args.tool](args)
