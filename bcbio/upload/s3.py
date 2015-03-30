"""Handle upload and retrieval of files from S3 on Amazon AWS.
"""
import datetime
import email
import os
import sys

from bcbio.distributed import objectstore
from bcbio.provenance import do
from bcbio.upload import filesystem

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

    region = "@%s" % config["region"] if config.get("region") else ""
    fname = "s3://%s%s/%s" % (config["bucket"], region, to_transfer[0][1])
    conn = objectstore.connect(fname)
    bucket = conn.lookup(config["bucket"])
    if not bucket:
        bucket = conn.create_bucket(config["bucket"])

    for fname, orig_keyname in to_transfer:
        keyname = os.path.join(config.get("folder", ""), orig_keyname)
        key = bucket.get_key(keyname) if bucket else None
        modified = datetime.datetime.fromtimestamp(email.utils.mktime_tz(
            email.utils.parsedate_tz(key.last_modified))) if key else None
        no_upload = key and modified >= finfo["mtime"]
        if not no_upload:
            if config.get("region") in objectstore.REGIONS_NEWPERMS["s3"]:
                _upload_file_aws_cli(fname, config["bucket"], keyname, config, finfo)
            else:
                _upload_file(fname, config["bucket"], keyname, config, finfo)

def _upload_file(fname, bucket, keyname, config=None, mditems=None):
    metadata = ["-m", "x-amz-server-side-encryption:AES256"]
    endpoint = []
    if mditems:
        for name, val in mditems.iteritems():
            val = _update_val(name, val)
            if val:
                metadata += ["-m", "x-amz-meta-%s:%s" % (name, val)]
    if config:
        if config.get("reduced_redundancy"):
            metadata += ["-m", "x-amz-storage-class:REDUCED_REDUNDANCY"]
        if config.get("region"):
            if config.get("region") != "us-east-1":
                endpoint = ["--endpoint=s3-%s.amazonaws.com" % config.get("region")]
    cmd = ["gof3r", "put", "--no-md5", "-b", bucket, "-k", keyname,
           "-p", fname] + endpoint + metadata
    do.run(cmd, "Upload to s3: %s %s" % (bucket, keyname))

def _upload_file_aws_cli(local_fname, bucket, keyname, config=None, mditems=None):
    """Potentially streaming download via the standard AWS command line interface.
    """
    s3_fname = "s3://%s/%s" % (bucket, keyname)
    args = ["--sse", "--expected-size", str(os.path.getsize(local_fname))]
    if config:
        if config.get("region"):
            args += ["--region", config.get("region")]
        if config.get("reduced_redundancy"):
            args += ["--storage-class", "REDUCED_REDUNDANCY"]
    cmd = [os.path.join(os.path.dirname(sys.executable), "aws"), "s3", "cp"] + args + \
          [local_fname, s3_fname]
    do.run(cmd, "Upload to s3: %s %s" % (bucket, keyname))

def upload_file_boto(fname, remote_fname, mditems=None):
    """Upload a file using boto instead of external tools.
    """
    r_fname = objectstore.parse_remote(remote_fname)
    conn = objectstore.connect(remote_fname)
    bucket = conn.lookup(r_fname.bucket)
    if not bucket:
        bucket = conn.create_bucket(r_fname.bucket)
    key = bucket.get_key(r_fname.key, validate=False)
    if mditems is None:
        mditems = {}
    if "x-amz-server-side-encryption" not in mditems:
        mditems["x-amz-server-side-encryption"] = "AES256"
    for name, val in mditems.iteritems():
        key.set_metadata(name, val)
    key.set_contents_from_filename(fname, encrypt_key=True)
