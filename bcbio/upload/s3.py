"""Handle upload and retrieval of files from S3 on Amazon AWS.
"""
import datetime
import email
import os

import boto

from bcbio.log import logger

def get_file(local_dir, bucket_name, fname, params):
    """Retrieve file from amazon S3 to a local directory for processing.
    """
    out_file = os.path.join(local_dir, os.path.basename(fname))
    conn = boto.connect_s3(params.get("access_key_id"), params.get("secret_access_key"))
    bucket = conn.get_bucket(bucket_name)
    key = bucket.get_key(fname)
    key.get_contents_to_filename(out_file)
    return out_file

def _update_val(key, val):
    if key == "mtime":
        return val.isoformat()
    elif key in ["path", "ext"]:
        return None
    else:
        return val

def update_file(finfo, sample_info, config):
    """Update the file to an Amazon S3 bucket.
    """
    conn = boto.connect_s3(config.get("access_key_id"),
                           config.get("secret_access_key"))
    bucket = conn.lookup(config["bucket"])
    if bucket is None:
        bucket = conn.create_bucket(config["bucket"])
    s3dirname = finfo["sample"] if finfo.has_key("sample") else finfo["run"]
    keyname = os.path.join(s3dirname, os.path.basename(finfo["path"]))
    key = bucket.get_key(keyname)
    modified = datetime.datetime.fromtimestamp(email.utils.mktime_tz(
        email.utils.parsedate_tz(key.last_modified))) if key else None
    no_upload = key and modified >= finfo["mtime"]
    if key is None:
        key = boto.s3.key.Key(bucket, keyname)
    if not no_upload:
        logger.info("Uploading to S3: %s %s" % (config["bucket"], keyname))
        for name, val in finfo.iteritems():
            val = _update_val(name, val)
            if val:
                key.set_metadata(name, val)
        key.set_contents_from_filename(finfo["path"],
          reduced_redundancy=config.get("reduced_redundancy", False))
