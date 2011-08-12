"""Top level management of analysis processing.

Handles copying remote files from sequencer, starting processing scripts,
and upload of results back to Galaxy.
"""
import os
import re
import contextlib

import yaml
import fabric.api as fabric
import fabric.contrib.files as fabric_files

from bcbio.pipeline import log
from bcbio.log import create_log_handler

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
        _, _, c_host_str = _config_hosts(config)
        with fabric.settings(host_string=c_host_str):
            fc_dir = _remote_copy(remote_info, config)
    return fc_dir

def _config_hosts(config):
    """Retrieve configured machines to perform analysis and copy on.
    """
    user = config["analysis"].get("user", None)
    host = config["analysis"].get("host", None)
    shell = config["analysis"].get("login_shell", None)
    if not user or not host:
        user = os.environ["USER"]
        host = re.sub(r'\..*', '', os.uname()[1])
    if not shell:
        shell = os.environ["SHELL"]
    copy_user = config["analysis"].get("copy_user", None)
    copy_host = config["analysis"].get("copy_host", None)
    if not copy_user or not copy_host:
        copy_user, copy_host = (user, host)
    analysis_host_str = "%s@%s" % (user, host)
    analysis_shell = "%s -i -l -c" % shell
    copy_host_str = "%s@%s" % (copy_user, copy_host)
    return analysis_host_str, analysis_shell, copy_host_str

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
    analysis_dir = os.path.join(config["analysis"]["base_dir"],
                                os.path.basename(remote_info["directory"]))
    with analysis_machine_config(config, config_file) as config_file:
        if not fabric_files.exists(analysis_dir):
            fabric.run("mkdir %s" % analysis_dir)
        with fabric.cd(analysis_dir):
            if config["algorithm"]["num_cores"] == "messaging":
                prog = config["analysis"].get("distributed_process_program",
                                              "distributed_nextgen_pipeline.py")
            else:
                prog = config["analysis"]["process_program"]
            cl = [prog, config_file, fc_dir]
            if run_yaml:
                cl.append(run_yaml)
            fabric.run(" ".join(cl))
    return analysis_dir

@contextlib.contextmanager
def analysis_machine_config(config, config_file):
    """Details about where to kick off the analysis and upload.
    """
    a_host_str, a_shell, _ = _config_hosts(config)
    if config["analysis"].get("config_file", None):
        config_file = config["analysis"]["config_file"]
    else:
        config_file = os.path.abspath(config_file)
    with fabric.settings(host_string=a_host_str, shell=a_shell):
        yield config_file

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
    with analysis_machine_config(config, config_file) as config_file:
        cl = [config["analysis"]["upload_program"], config_file, fc_dir,
              analysis_dir]
        if run_yaml:
            cl.append(run_yaml)
        fabric.run(" ".join(cl))
