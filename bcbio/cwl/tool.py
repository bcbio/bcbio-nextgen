"""Run bcbio generated CWL with a supported tool.

Handles wrapping and integrating with multiple tools making it easier
to run bcbio in a standard way in many environments.
"""
from __future__ import print_function
import glob
import json
import os
import shutil
import subprocess
import sys

from bcbio import utils
from bcbio.cwl import hpc
from bcbio.distributed import objectstore

def _get_main_and_json(directory):
    """Retrieve the main CWL and sample JSON files from a bcbio generated directory.
    """
    directory = os.path.normpath(os.path.abspath(directory))
    checker_main = os.path.normpath(os.path.join(directory, os.path.pardir, "checker-workflow-wrapping-tool.cwl"))
    if checker_main and os.path.exists(checker_main):
        main_cwl = [checker_main]
    else:
        main_cwl = glob.glob(os.path.join(directory, "main-*.cwl"))
        main_cwl = [x for x in main_cwl if not x.find("-pack") >= 0]
        assert len(main_cwl) == 1, "Did not find main CWL in %s" % directory
    main_json = glob.glob(os.path.join(directory, "main-*-samples.json"))
    assert len(main_json) == 1, "Did not find main json in %s" % directory
    project_name = os.path.basename(directory).split("-workflow")[0]
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

def _pack_cwl(unpacked_cwl):
    """Pack CWL into a single document for submission.
    """
    out_file = "%s-pack%s" % os.path.splitext(unpacked_cwl)
    cmd = "cwltool --pack {unpacked_cwl} > {out_file}"
    _run_tool(cmd.format(**locals()))
    return out_file

def _chown_workdir(work_dir):
    """Ensure work directory files owned by original user.

    Docker runs can leave root owned files making cleanup difficult.
    Skips this if it fails, avoiding errors where we run remotely
    and don't have docker locally.
    """
    cmd = ("""docker run --rm -v %s:%s quay.io/bcbio/bcbio-base /bin/bash -c 'chown -R %s %s'""" %
           (work_dir, work_dir, os.getuid(), work_dir))
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError:
        pass

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

def _run_wes(args):
    """Run CWL using a Workflow Execution Service (WES) endpoint
    """
    main_file, json_file, project_name = _get_main_and_json(args.directory)
    main_file = _pack_cwl(main_file)
    if args.host and "stratus" in args.host:
        _run_wes_stratus(args, main_file, json_file)
    else:
        opts = ["--no-wait"]
        if args.host:
            opts += ["--host", args.host]
        if args.auth:
            opts += ["--auth", args.auth]
        cmd = ["wes-client"] + opts + [main_file, json_file]
        _run_tool(cmd)

def _run_wes_stratus(args, main_file, json_file):
    """Run WES on Illumina stratus endpoint server, which wes-client doesn't support.

    https://stratus-docs.readme.io/docs/quick-start-4
    """
    import requests
    base_url = args.host
    if not base_url.startswith("http"):
        base_url = "https://%s" % base_url
    with open(main_file) as in_handle:
        r = requests.post("%s/v1/workflows" % base_url,
                          headers={"Content-Type": "application/json",
                                   "Authorization": "Bearer %s" % args.auth},
                          data=in_handle.read())
    print(r.status_code)
    print(r.text)

def _estimate_runner_memory(json_file):
    """Estimate Java memory requirements based on number of samples.

    A rough approach to selecting correct allocated memory for Cromwell.
    """
    with open(json_file) as in_handle:
        sinfo = json.load(in_handle)
    num_parallel = 1
    for key in ["config__algorithm__variantcaller", "description"]:
        item_counts = []
        n = 0
        for val in (sinfo.get(key) or []):
            n += 1
            if val:
                if isinstance(val, (list, tuple)):
                    item_counts.append(len(val))
                else:
                    item_counts.append(1)
        print(key, n, item_counts)
        if n and item_counts:
            num_parallel = n * max(item_counts)
            break
    if num_parallel < 25:
        return "3g"
    if num_parallel < 150:
        return "6g"
    elif num_parallel < 500:
        return "12g"
    else:
        return "24g"

def _run_cromwell(args):
    """Run CWL with Cromwell.
    """
    main_file, json_file, project_name = _get_main_and_json(args.directory)
    work_dir = utils.safe_makedir(os.path.join(os.getcwd(), "cromwell_work"))
    final_dir = utils.safe_makedir(os.path.join(work_dir, "final"))
    if args.no_container:
        _remove_bcbiovm_path()
    log_file = os.path.join(work_dir, "%s-cromwell.log" % project_name)
    metadata_file = os.path.join(work_dir, "%s-metadata.json" % project_name)
    option_file = os.path.join(work_dir, "%s-options.json" % project_name)
    cromwell_opts = {"final_workflow_outputs_dir": final_dir,
                     "default_runtime_attributes": {"bootDiskSizeGb": 20}}
    with open(option_file, "w") as out_handle:
        json.dump(cromwell_opts, out_handle)

    cmd = ["cromwell", "-Xms1g", "-Xmx%s" % _estimate_runner_memory(json_file),
           "run", "--type", "CWL",
           "-Dconfig.file=%s" % hpc.create_cromwell_config(args, work_dir, json_file)]
    cmd += hpc.args_to_cromwell_cl(args)
    cmd += ["--metadata-output", metadata_file, "--options", option_file,
            "--inputs", json_file, main_file]
    with utils.chdir(work_dir):
        _run_tool(cmd, not args.no_container, work_dir, log_file)
        if metadata_file and utils.file_exists(metadata_file):
            with open(metadata_file) as in_handle:
                metadata = json.load(in_handle)
            if metadata["status"] == "Failed":
                _cromwell_debug(metadata)
                sys.exit(1)
            else:
                _cromwell_move_outputs(metadata, final_dir)

def _cromwell_debug(metadata):
    """Format Cromwell failures to make debugging easier.
    """
    def get_failed_calls(cur, key=None):
        if key is None: key = []
        out = []
        if isinstance(cur, dict) and "failures" in cur and "callRoot" in cur:
            out.append((key, cur))
        elif isinstance(cur, dict):
            for k, v in cur.items():
                out.extend(get_failed_calls(v, key + [k]))
        elif isinstance(cur, (list, tuple)):
            for i, v in enumerate(cur):
                out.extend(get_failed_calls(v, key + [i]))
        return out
    print("Failed bcbio Cromwell run")
    print("-------------------------")
    for fail_k, fail_call in get_failed_calls(metadata["calls"]):
        root_dir = os.path.join("cromwell_work", os.path.relpath(fail_call["callRoot"]))
        print("Failure in step: %s" % ".".join([str(x) for x in fail_k]))
        print("  bcbio log file     : %s" % os.path.join(root_dir, "execution", "log", "bcbio-nextgen-debug.log"))
        print("  bcbio commands file: %s" % os.path.join(root_dir, "execution", "log",
                                                         "bcbio-nextgen-commands.log"))
        print("  Cromwell directory : %s" % root_dir)
        print()

def _cromwell_move_outputs(metadata, final_dir):
    """Move Cromwell outputs to the final upload directory.
    """
    sample_key = [k for k in metadata["outputs"].keys() if k.endswith(("rgnames__sample", "rgnames__sample_out"))][0]
    project_dir = utils.safe_makedir(os.path.join(final_dir, "project"))
    samples = metadata["outputs"][sample_key]
    def _copy_with_secondary(f, dirname):
        if len(f["secondaryFiles"]) > 1:
            dirname = utils.safe_makedir(os.path.join(dirname, os.path.basename(os.path.dirname(f["location"]))))
        if not objectstore.is_remote(f["location"]):
            finalf = os.path.join(dirname, os.path.basename(f["location"]))
            if not utils.file_uptodate(finalf, f["location"]):
                shutil.copy(f["location"], dirname)
        [_copy_with_secondary(sf, dirname) for sf in f["secondaryFiles"]]
    def _write_to_dir(val, dirname):
        if isinstance(val, (list, tuple)):
            [_write_to_dir(v, dirname) for v in val]
        else:
            _copy_with_secondary(val, dirname)
    for k, vals in metadata["outputs"].items():
        if k != sample_key:
            if k.endswith(("summary__multiqc")):
                vs = [v for v in vals if v]
                assert len(vs) == 1
                _write_to_dir(vs[0], project_dir)
            elif len(vals) == len(samples):
                for s, v in zip(samples, vals):
                    if v:
                        _write_to_dir(v, utils.safe_makedir(os.path.join(final_dir, s)))
            elif len(vals) == 1:
                _write_to_dir(vals[0], project_dir)
            elif len(vals) > 0:
                raise ValueError("Unexpected sample and outputs: %s %s %s" % (k, samples, vals))

def _run_sbgenomics(args):
    """Run CWL on SevenBridges platform and Cancer Genomics Cloud.
    """
    assert not args.no_container, "Seven Bridges runs require containers"
    main_file, json_file, project_name = _get_main_and_json(args.directory)
    flags = []
    cmd = ["sbg-cwl-runner"] + flags + args.toolargs + [main_file, json_file]
    _run_tool(cmd)

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
          "cromwell": _run_cromwell,
          "arvados": _run_arvados,
          "sbg": _run_sbgenomics,
          "toil": _run_toil,
          "bunny": _run_bunny,
          "funnel": _run_funnel,
          "wes": _run_wes}

def run(args):
    _TOOLS[args.tool](args)
