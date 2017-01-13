"""
    Handle upload and retrieval of files from iRODS. This method requires a preconfigured connection
    to an iRODS repository throught the `iinit` command
    config options:
    upload:
      method: irods
      dir: absolute parent path in iRODS repository
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
<<<<<<< HEAD
    meta_string=[]
    #assuming metadata is a dictionary: {key1:value1,key2:value2}
    for keyname, value in metadata.iteritems():
        unit=""
        meta_string += [str(keyname).replace(" ","_"),str(value).replace(" ","_"),unit]
    return ";".join(meta_string)
=======
    meta_string=""
    #assuming metadata is a dictionary: {key1:value1,key2:value2}
    for keyname, value in metadata.iteritems():
        meta_string += str(keyname;value;;)
    return meta_string
>>>>>>> d94647cc433bef6b72ba13fdc6ce58ff79a654d0


def update_file(finfo, sample_info, config):
    """
    Update the file to an iRODS repository.
    """
    ffinal = filesystem.update_file(finfo, sample_info, config, pass_uptodate=True)
    if os.path.isdir(ffinal):
        to_transfer = []
        for path, dirs, files in os.walk(ffinal):
            for f in files:
                full_f = os.path.join(path, f)
                k = full_f.replace(os.path.abspath(config["dir"]) + "/", "")
                to_transfer.append((full_f, k))
    else:
        k = ffinal.replace(os.path.abspath(config["dir"]) + "/", "")
        to_transfer = [(ffinal, k)]
    fname = "%s/%s" % (config["dir"], to_transfer[0][1])

    for fname, orig_keyname in to_transfer:
        keyname = os.path.join(config.get("dir", ""), orig_keyname)
<<<<<<< HEAD
        #metadata= _format_metadata(sample_info['summary']['metrics'])
        _upload_file_icommands_cli(fname, keyname, config )
=======
        metadata= _format_metadata(sample_info)
        if not no_upload:
            _upload_file_icommands_cli(fname, keyname, config, metadata)
>>>>>>> d94647cc433bef6b72ba13fdc6ce58ff79a654d0

def _upload_file_icommands_cli(local_fname, keyname, config=None, metadata=None):
    """
    Upload via the standard icommands CLI.
    example: iput -K -P -R $resource -v --metadata "attr1;val1;unit1;attr2;val2;unit2;" $file $dir/$file
    go to https://docs.irods.org/4.2.0/icommands/user/#iput for more info
    """

    irods_fname = "%s/%s" % (config.get("dir"), keyname)
    args = ["-K","-P", "-v"]
    if config:
        if config.get("resource"):
            args += ["-R", config.get("resource")]
        if config.get("ticket"):
            args += ["-t", config.get("ticket")]
        if config.get("extra"):
            args += config.get("extra")
    if metadata:
        args += ["--metadata", metadata]

    cmd = ["iput"] + args + [local_fname, irods_fname]
    do.run(cmd, "Upload to iRODS: %s" % (keyname))

def download_file(file):
    """
    Download file from an iRODS repository
    """
    _download_file_icommands_cli(file)

def _download_file_icommands_cli(file):
    """
    Download file from iRODS using the iCommands CLI
    example: iget -K $file
    """
    cmd = ["iget","-K"] + [file]
    do.run(cmd, "Download from iRODS: %s" % (file))
