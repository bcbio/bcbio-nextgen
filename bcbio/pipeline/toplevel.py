"""Top level management of analysis processing.

Handles copying remote files from sequencer, starting processing scripts,
and upload of results back to Galaxy.
"""
import os
import re
import subprocess

import yaml
# Fabric only needed on running side, not on setup and initial import
try:
    import fabric.api as fabric
    import fabric.contrib.files as fabric_files
except ImportError:
    fabric, fabric_files = (None, None)

from bcbio.pipeline import log
from bcbio.log import create_log_handler
from bcbio import utils

def analyze_and_upload(remote_info, config_file):
    """Main entry point for analysis and upload to Galaxy.
    """
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    log_handler = create_log_handler(config, log.name)
    with log_handler.applicationbound():
        fc_dir = _copy_from_sequencer(remote_info, config)
        analysis_dir = _run_analysis(fc_dir, remote_info, config, config_file)
        _upload_to_galaxy(fc_dir, analysis_dir, remote_info,
                          config, config_file)

# ## Copying over files from sequencer, if necessary

def _copy_from_sequencer(remote_info, config):
    """Get local directory of flowcell info, or copy from sequencer.
    """
    if remote_info.has_key("fc_dir"):
        fc_dir = remote_info["fc_dir"]
        assert os.path.exists(fc_dir)
    else:
        log.debug("Remote host information: %s" % remote_info)
        c_host_str = _config_hosts(config)
        with fabric.settings(host_string=c_host_str):
            fc_dir = _remote_copy(remote_info, config)
    return fc_dir

def _config_hosts(config):
    """Retrieve configured machines to perform analysis and copy on.
    """
    copy_user = config["analysis"].get("copy_user", None)
    copy_host = config["analysis"].get("copy_host", None)
    if not copy_user or not copy_host:
        copy_user = os.environ["USER"]
        copy_host = re.sub(r'\..*', '', os.uname()[1])
    copy_host_str = "%s@%s" % (copy_user, copy_host)
    return copy_host_str

def _remote_copy(remote_info, config):
    """Securely copy files from remote directory to the processing server.

    This requires ssh public keys to be setup so that no password entry
    is necessary.
    """
    fc_dir = os.path.join(config["analysis"]["store_dir"],
                          os.path.basename(remote_info['directory']))
    log.info("Copying analysis files to %s" % fc_dir)
    if not fabric_files.exists(fc_dir):
        fabric.run("mkdir %s" % fc_dir)
    for fcopy in remote_info['to_copy']:
        target_loc = os.path.join(fc_dir, fcopy)
        if not fabric_files.exists(target_loc):
            target_dir = os.path.dirname(target_loc)
            if not fabric_files.exists(target_dir):
                fabric.run("mkdir -p %s" % target_dir)
            cl = ["scp", "-r", "%s@%s:%s/%s" %
                  (remote_info["user"], remote_info["hostname"],
                   remote_info["directory"], fcopy),
                  target_loc]
            fabric.run(" ".join(cl))
    log.info("Analysis files copied")
    return fc_dir

def _run_analysis(fc_dir, remote_info, config, config_file):
    """Run local or distributed analysis, wait to finish.
    """
    run_yaml = _get_run_yaml(remote_info, fc_dir, config)
    analysis_dir = os.path.join(config["analysis"].get("base_dir", os.getcwd()),
                                os.path.basename(remote_info["directory"]))
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
    with utils.chdir(analysis_dir):
        if config["algorithm"]["num_cores"] == "messaging":
            prog = config["analysis"].get("distributed_process_program",
                                          "distributed_nextgen_pipeline.py")
        else:
            prog = config["analysis"]["process_program"]
        cl = [prog, config_file, fc_dir]
        if run_yaml:
            cl.append(run_yaml)
        subprocess.check_call(cl)
    return analysis_dir

def _get_run_yaml(remote_info, fc_dir, config):
    """Retrieve YAML specifying run from configured or default location.
    """
    if remote_info.get("run_yaml", None):
        run_yaml = remote_info["run_yaml"]
    else:
        run_yaml = os.path.join(config["analysis"]["store_dir"],
                                os.path.basename(fc_dir), "run_info.yaml")
    if not os.path.exists(run_yaml):
        run_yaml = None
    return run_yaml

def _upload_to_galaxy(fc_dir, analysis_dir, remote_info, config, config_file):
    """Upload results from analysis directory to Galaxy data libraries.
    """
    run_yaml = _get_run_yaml(remote_info, fc_dir, config)
    with utils.chdir(analysis_dir):
        cl = [config["analysis"]["upload_program"], config_file, fc_dir,
              analysis_dir]
        if run_yaml:
            cl.append(run_yaml)
        subprocess.check_call(cl)
