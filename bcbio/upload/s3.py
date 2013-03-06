"""Handle upload and retrieval of files from S3 on Amazon AWS.
"""
import os

import boto

def get_file(local_dir, bucket_name, fname, params):
    """Retrieve file from amazon S3 to a local directory for processing.
    """
    out_file = os.path.join(local_dir, os.path.basename(fname))
    conn = boto.connect_s3(params.get("access_key_id"), params.get("secret_access_key"))
    bucket = conn.get_bucket(bucket_name)
    key = bucket.get_key(fname)
    key.get_contents_to_filename(out_file)
    return out_file

def update_file(finfo, sample_info, config):
    """Update the file to an Amazon S3 bucket.
    """
    raise NotImplementedError("TODO")
