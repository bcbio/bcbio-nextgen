"""
    Handle upload and retrieval of files from iRODS. This method requires a preconfigured connection
    to an iRODS repository throught the `iinit` command
    config options:
    upload:
      dir: ../final
      method: irods
      folder: absolute parent path in iRODS repository
      resource: (optional) iRODS resource name, if other than default
"""
import os

from bcbio.provenance import do
from bcbio.upload import filesystem

def _check_create_collection(irods_fname,isdir=False):
        irods_dir=""
        if isdir:
            irods_dir=irods_fname
        else:
            irods_dir=os.path.dirname(irods_fname)
        cmd = ["imkdir", "-p",irods_dir]
        do.run(cmd,"iRODS: create collection %s" % (irods_dir))

def update_file(finfo, sample_info, config):
    """
    Update the file to an iRODS repository.
    """
    ffinal = filesystem.update_file(finfo, sample_info, config, pass_uptodate=True)

    _upload_dir_icommands_cli(config.get("dir"), config.get("folder"), config)

def _upload_dir_icommands_cli(local_dir, irods_dir, config=None, metadata=None):
    """
    Upload directory recursively via the standard icommands CLI.
    example: irsync -Kvar -R $resource $local_dir i:$irods_dir
    go to https://docs.irods.org/4.2.0/icommands/user/#irsync for more info
    """

    args = ["-K","-v","-a","-r"]
    if config:
        if config.get("resource"):
            args += ["-R", config.get("resource")]

    _check_create_collection(irods_dir,isdir=True)
    cmd = ["irsync"] + args + [local_dir, "i:"+irods_dir]
    do.run(cmd, "Uploading to iRODS")


