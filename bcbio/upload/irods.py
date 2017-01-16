"""
    Handle upload and retrieval of files from iRODS. This method requires a preconfigured connection
    to an iRODS repository throught the `iinit` command
    config options:
    upload:
      dir: ../final
      method: irods
      folder: absolute parent path in iRODS repository
      resource: (optional) iRODS resource name, if other than default
      ticket: (optional) iRODS ticket, for ticket based access
      extra: (optional)["list","of","arbitrary","options","that","can","be","passed","to","iput"]

"""
import datetime
import email
import os
import sys

from bcbio.provenance import do
from bcbio.upload import filesystem

def _format_metadata(metadata):
    """
    Format metadata to use in iCommands CLI
    requisite format: "attr1;val1;unit1;attr2;val2;unit2;"
    """
    meta_string=[]
    #assuming metadata is a dictionary: {key1:value1,key2:value2}
    for keyname, value in metadata.iteritems():
        unit=""
        meta_string += [str(keyname).replace(" ","_"),str(value).replace(" ","_"),unit]
    return ";".join(meta_string)

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

    #Loop over files and upload file by file - ineffecient - lots of subshells
    #if os.path.isdir(ffinal):
    #    to_transfer = []
    #    for path, dirs, files in os.walk(ffinal):
    #        for f in files:
    #            full_f = os.path.join(path, f)
    #            k = full_f.replace(os.path.abspath(config["dir"]) + "/", "")
    #            to_transfer.append((full_f, k))
    #else:
    #    k = ffinal.replace(os.path.abspath(config["dir"]) + "/", "")
    #    to_transfer = [(ffinal, k)]
    #for fname, orig_keyname in to_transfer:
    #    keyname = os.path.join(config.get("folder", ""), orig_keyname)
    #    metadata= _format_metadata(sample_info['summary']['metrics'])
    #    _upload_file_icommands_cli(fname, keyname, config )

def _upload_dir_icommands_cli(local_dir, irods_dir, config=None, metadata=None):
    """
    Upload directory recursively via the standard icommands CLI.
    example: iput -KPvr -R $resource --metadata "attr1;val1;unit1;attr2;val2;unit2;" $file $dir/$file
    go to https://docs.irods.org/4.2.0/icommands/user/#iput for more info
    """

    args= ["-K","-T","-v","-r","-f","-X","%s/upload.iRODS"%(config.get("dir")),"--lfrestart","%s/lfrestart.iRODS"%(config.get("dir")),"--retries 5"]
    if config:
        if config.get("resource"):
            args += ["-R", config.get("resource")]
        if config.get("ticket"):
            args += ["-t", config.get("ticket")]
        if config.get("extra"):
            args += config.get("extra")
    if metadata:
        args += ["--metadata", metadata]

    _check_create_collection(irods_dir,isdir=True)
    cmd = ["iput"] + args + [local_dir.rstrip("/"), irods_dir.rstrip("/")]
    do.run(cmd, "Uploading to iRODS")

def _upload_file_icommands_cli(local_fname, irods_fname, config=None, metadata=None):
    """
    Upload via the standard icommands CLI.
    example: iput -K -P -R $resource -v --metadata "attr1;val1;unit1;attr2;val2;unit2;" $file $dir/$file
    go to https://docs.irods.org/4.2.0/icommands/user/#iput for more info
    """

    args = ["-K","-T","-v","-f","-X","%s/upload.iRODS"%(config.get("dir")),"--lfrestart","%s/lfrestart.iRODS"%(config.get("dir")),"--retries 5"]
    if config:
        if config.get("resource"):
            args += ["-R", config.get("resource")]
        if config.get("ticket"):
            args += ["-t", config.get("ticket")]
        if config.get("extra"):
            args += config.get("extra")
    if metadata:
        args += ["--metadata", metadata]

    _check_create_collection(irods_fname)

    cmd = ["iput"] + args + [local_fname, irods_fname]
    do.run(cmd, "Upload to iRODS: %s" % (irods_fname))
