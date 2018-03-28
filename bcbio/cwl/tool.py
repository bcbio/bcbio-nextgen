"""Run bcbio generated CWL with a supported tool.

Handles wrapping and integrating with multiple tools making it easier
to run bcbio in a standard way in many environments.
"""
import glob
import os
import shutil
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
    project_name = os.path.basename(directory).replace("-workflow", "")
    return main_cwl[0], main_json[0], project_name

def _run_tool(cmd, use_container=True, work_dir=None, log_file=None):
    """Run with injection of bcbio path.

    Place at end for runs without containers to avoid overriding other
    bcbio installations.
    """
    if isinstance(cmd, (list, tuple)):
        cmd = " ".join([str(x) for x in cmd])
    cmd = utils.local_path_export(at_start=use_container) + cmd
    if log_file:
        cmd += " 2>&1 | tee -a %s" % log_file
    try:
        print("Running: %s" % cmd)
        subprocess.check_call(cmd, shell=True)
    finally:
        if use_container and work_dir:
            _chown_workdir(work_dir)

def _chown_workdir(work_dir):
    """Ensure work directory files owned by original user.

    Docker runs can leave root owned files making cleanup difficult.
    """
    cmd = ("""docker run --rm -v %s:%s quay.io/bcbio/bcbio-base /bin/bash -c 'chown -R %s %s'""" %
           (work_dir, work_dir, os.getuid(), work_dir))
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
    main_file, json_file, project_name = _get_main_and_json(args.directory)
    work_dir = utils.safe_makedir(os.path.join(os.getcwd(), "cwltool_work"))
    tmp_dir = utils.safe_makedir(os.path.join(work_dir, "tmpcwl"))
    log_file = os.path.join(work_dir, "%s-cwltool.log" % project_name)
    os.environ["TMPDIR"] = tmp_dir
    flags = ["--tmpdir-prefix", tmp_dir, "--tmp-outdir-prefix", tmp_dir]
    if args.no_container:
        _remove_bcbiovm_path()
        flags += ["--no-container", "--preserve-environment", "PATH", "--preserve-environment", "HOME"]
    cmd = ["cwltool"] + flags + args.toolargs + ["--", main_file, json_file]
    with utils.chdir(work_dir):
        _run_tool(cmd, not args.no_container, work_dir, log_file=log_file)

def _run_arvados(args):
    """Run CWL on Arvados.
    """
    assert not args.no_container, "Arvados runs require containers"
    assert "ARVADOS_API_TOKEN" in os.environ and "ARVADOS_API_HOST" in os.environ, \
        "Need to set ARVADOS_API_TOKEN and ARVADOS_API_HOST in environment to run"
    main_file, json_file, project_name = _get_main_and_json(args.directory)
    flags = ["--enable-reuse", "--api", "containers", "--submit", "--no-wait"]
    cmd = ["arvados-cwl-runner"] + flags + args.toolargs + [main_file, json_file]
    _run_tool(cmd)

def _run_toil(args):
    """Run CWL with Toil.
    """
    main_file, json_file, project_name = _get_main_and_json(args.directory)
    work_dir = utils.safe_makedir(os.path.join(os.getcwd(), "toil_work"))
    tmp_dir = utils.safe_makedir(os.path.join(work_dir, "tmpdir"))
    os.environ["TMPDIR"] = tmp_dir
    log_file = os.path.join(work_dir, "%s-toil.log" % project_name)
    jobstore = os.path.join(work_dir, "cwltoil_jobstore")
    flags = ["--jobStore", jobstore, "--logFile", log_file, "--workDir", tmp_dir, "--linkImports"]
    if os.path.exists(jobstore):
        flags += ["--restart"]
    # caching causes issues for batch systems
    if "--batchSystem" in args.toolargs:
        flags += ["--disableCaching"]
    flags += args.toolargs
    if args.no_container:
        _remove_bcbiovm_path()
        flags += ["--no-container", "--preserve-environment", "PATH", "HOME"]
    cmd = ["cwltoil"] + flags + ["--", main_file, json_file]
    with utils.chdir(work_dir):
        _run_tool(cmd, not args.no_container, work_dir)
        for tmpdir in (glob.glob(os.path.join(work_dir, "out_tmpdir*")) +
                       glob.glob(os.path.join(work_dir, "tmp*"))):
            if os.path.isdir(tmpdir):
                shutil.rmtree(tmpdir)

def _run_bunny(args):
    """Run CWL with rabix bunny.
    """
    main_file, json_file, project_name = _get_main_and_json(args.directory)
    work_dir = utils.safe_makedir(os.path.join(os.getcwd(), "bunny_work"))
    flags = ["-b", work_dir]
    log_file = os.path.join(work_dir, "%s-bunny.log" % project_name)
    if os.path.exists(work_dir):
        caches = [os.path.join(work_dir, d) for d in os.listdir(work_dir)
                  if os.path.isdir(os.path.join(work_dir, d))]
        if caches:
            flags += ["--cache-dir", max(caches, key=os.path.getmtime)]
    if args.no_container:
        _remove_bcbiovm_path()
        flags += ["--no-container"]
    cmd = ["rabix"] + flags + [main_file, json_file]
    with utils.chdir(work_dir):
        _run_tool(cmd, not args.no_container, work_dir, log_file)

def _run_funnel(args):
    """Run funnel TES server with rabix bunny for CWL.
    """
    host = "localhost"
    port = "8088"
    main_file, json_file, project_name = _get_main_and_json(args.directory)
    work_dir = utils.safe_makedir(os.path.join(os.getcwd(), "funnel_work"))
    log_file = os.path.join(work_dir, "%s-funnel.log" % project_name)
    # Create bunny configuration directory with TES backend
    orig_config_dir = os.path.join(os.path.dirname(os.path.realpath(utils.which("rabix"))), "config")
    work_config_dir = utils.safe_makedir(os.path.join(work_dir, "rabix_config"))
    for fname in os.listdir(orig_config_dir):
        if fname == "core.properties":
            with open(os.path.join(orig_config_dir, fname)) as in_handle:
                with open(os.path.join(work_config_dir, fname), "w") as out_handle:
                    for line in in_handle:
                        if line.startswith("backend.embedded.types"):
                            line = "backend.embedded.types=TES\n"
                        out_handle.write(line)
        else:
            shutil.copy(os.path.join(orig_config_dir, fname), os.path.join(work_config_dir, fname))
    flags = ["-c", work_config_dir,
             "-tes-url=http://%s:%s" % (host, port), "-tes-storage=%s" % work_dir]
    if args.no_container:
        _remove_bcbiovm_path()
        flags += ["--no-container"]
    cmd = ["rabix"] + flags + [main_file, json_file]
    funnelp = subprocess.Popen(["funnel", "server", "run",
                                "--Server.HostName", host, "--Server.HTTPPort", port,
                                "--LocalStorage.AllowedDirs", work_dir,
                                "--Worker.WorkDir", os.path.join(work_dir, "funnel-work")])
    try:
        with utils.chdir(work_dir):
            _run_tool(cmd, not args.no_container, work_dir, log_file)
    finally:
        funnelp.kill()

_TOOLS = {"cwltool": _run_cwltool,
          "arvados": _run_arvados,
          "toil": _run_toil,
          "bunny": _run_bunny,
          "funnel": _run_funnel}

def run(args):
    _TOOLS[args.tool](args)
