"""Handle upload and retrieval of files from S3 on Amazon AWS.
"""
import datetime
import email
import os

import boto

from bcbio import utils
from bcbio.provenance import do

def get_file(local_dir, bucket_name, fname, params):
    """Retrieve file from amazon S3 to a local directory for processing.
    """
    out_file = os.path.join(local_dir, os.path.basename(fname))
    if not utils.file_exists(out_file):
        metadata = []
        if params.get("reduced_redundancy"):
            metadata += ["-m", "x-amz-storage-class:REDUCED_REDUNDANCY"]
        cmd = ["gof3r", "get", "--no-md5", "-b", bucket_name, "-k", fname,
               "-p", out_file] + metadata
        do.run(cmd, "Retrieve from s3")
    return out_file

def _update_val(key, val):
    if key == "mtime":
        return val.isoformat()
    elif key in ["path", "ext"]:
        return None
    else:
        return val

def update_file(finfo, sample_info, config):
    """Update the file to an Amazon S3 bucket, using server side encryption.
    """
    conn = boto.connect_s3()
    s3dirname = finfo["sample"] if "sample" in finfo else finfo["run"]
    keyname = os.path.join(s3dirname, os.path.basename(finfo["path"]))

    bucket = conn.lookup(config["bucket"])
    key = bucket.get_key(keyname) if bucket else None
    modified = datetime.datetime.fromtimestamp(email.utils.mktime_tz(
        email.utils.parsedate_tz(key.last_modified))) if key else None
    no_upload = key and modified >= finfo["mtime"]
    if not no_upload:
        metadata = ["-m", "x-amz-server-side-encryption:AES256"]
        for name, val in finfo.iteritems():
            val = _update_val(name, val)
            if val:
                metadata += ["-m", "-m x-amz-meta-%s:%s" % (key, val)]
        cmd = ["gof3r", "put", "--no-md5", "-b", config["bucket"], "-k", keyname,
               "-p", finfo["path"]] + metadata
        do.run(cmd, "Upload to s3: %s %s" % (config["bucket"], keyname))
