"""Transfer raw files from finished NGS runs for backup and storage.
"""
import os

import yaml
import fabric.api as fabric
import fabric.contrib.files as fabric_files

from bcbio.pipeline import log
from bcbio.log import create_log_handler

def long_term_storage(remote_info, config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    log_handler = create_log_handler(config, log.name)
    with log_handler.applicationbound():
        _copy_for_storage(remote_info, config)

def _copy_for_storage(remote_info, config):
    """Securely copy files from remote directory to the storage server.

    This requires ssh public keys to be setup so that no password entry
    is necessary, Fabric is used to manage setting up copies on the remote
    storage server.
    """
    log.info("Copying run data over to remote storage: %s" % config["store_host"])
    log.debug("The contents from AMQP for this dataset are:\n %s" % remote_info)
    base_dir = config["store_dir"]
    fabric.env.host_string = "%s@%s" % (config["store_user"], config["store_host"])
    fc_dir = os.path.join(base_dir, os.path.basename(remote_info['directory']))
    if not fabric_files.exists(fc_dir):
        fabric.run("mkdir %s" % fc_dir)
    for fcopy in remote_info['to_copy']:
        target_loc = os.path.join(fc_dir, fcopy)
        if not fabric_files.exists(target_loc):
            target_dir = os.path.dirname(target_loc)
            if not fabric_files.exists(target_dir):
                fabric.run("mkdir -p %s" % target_dir)
            cl = ["scp", "-r", "%s@%s:%s/%s" % (
                  remote_info["user"], remote_info["hostname"], remote_info["directory"],
                  fcopy), target_loc]
            fabric.run(" ".join(cl))
