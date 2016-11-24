"""Run bcbio generated CWL with a supported tool.

Handles wrapping and integrating with multiple tools making it easier
to run bcbio in a standard way in many environments.
"""
import glob
import os
import subprocess
import sys

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

def _run_tool(cmd, use_container=True):
    """Run with injection of bcbio path.

    Place at end for runs without containers to avoid overriding other
    bcbio installations.
    """
    if isinstance(cmd, (list, tuple)):
        cmd = " ".join([str(x) for x in cmd])
    cmd = utils.local_path_export(at_start=use_container) + cmd
    subprocess.check_call(cmd, shell=True)

def _remove_bcbiovm_path():
    """Avoid referencing minimal bcbio_nextgen in bcbio_vm installation.
    """
    cur_path = os.path.dirname(os.path.realpath(sys.executable))
    paths = os.environ["PATH"].split(":")
    if cur_path in paths:
        paths.remove(cur_path)
        os.environ["PATH"] = ":".join(paths)

def _run_cwltool(args):
    """Run with cwltool -- reference implementation.
    """
    main_file, json_file = _get_main_and_json(args.directory)
    work_dir = utils.safe_makedir(os.path.join(os.getcwd(), "cwltool_work"))
    tmp_dir = utils.safe_makedir(os.path.join(work_dir, "tmpcwl"))
    os.environ["TMPDIR"] = tmp_dir
    flags = ["--tmpdir-prefix", tmp_dir, "--tmp-outdir-prefix", tmp_dir]
    if args.no_container:
        _remove_bcbiovm_path()
        flags += ["--no-container", "--preserve-environment", "PATH", "--preserve-environment", "HOME"]
    cmd = ["cwltool"] + flags + args.toolargs + ["--", main_file, json_file]
    _run_tool(cmd, not args.no_container)

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
    os.environ["TMPDIR"] = work_dir
    log_file = os.path.join(work_dir, "cwltoil.log")
    jobstore = os.path.join(work_dir, "cwltoil_jobstore")
    flags = ["--jobStore", jobstore, "--logFile", log_file, "--workDir", work_dir]
    if os.path.exists(jobstore):
        flags += ["--restart"]
    if "--batchSystem" in args.toolargs:
        flags += ["--disableCaching"]
    flags += args.toolargs
    if args.no_container:
        _remove_bcbiovm_path()
        flags += ["--no-container", "--preserve-environment", "PATH", "HOME"]
    cmd = ["cwltoil"] + flags + ["--", main_file, json_file]
    with utils.chdir(work_dir):
        _run_tool(cmd, not args.no_container)

_TOOLS = {"cwltool": _run_cwltool,
          "arvados": _run_arvados,
          "toil": _run_toil}

def run(args):
    _TOOLS[args.tool](args)
