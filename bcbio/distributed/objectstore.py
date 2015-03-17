"""Manage pushing and pulling files from an object store like Amazon Web Services S3.
"""
import collections
import os
import subprocess
import zlib

import boto

from bcbio.distributed.transaction import file_transaction
from bcbio import utils

# ## definitions

SUPPORTED_REMOTES = ("s3://",)
BIODATA_INFO = {"S3": "s3://biodata/prepped/{build}/{build}-{target}.tar.gz"}

# ## Utilities

def is_remote(fname):
    return isinstance(fname, basestring) and fname.lower().startswith(SUPPORTED_REMOTES)

def file_exists_or_remote(fname):
    """Check if a file exists or is accessible remotely.
    """
    if is_remote(fname):
        return True
    else:
        return utils.file_exists(fname)

def parse_remote(fname):
    """Parses a remote filename into bucket and key information.

    Handles S3 with optional region name specified in key:
      BUCKETNAME@REGIONNAME/KEY
    """
    RemoteFile = collections.namedtuple("RemoteFile", "bucket,key,region")
    if fname.startswith("s3://"):
        bucket, key = fname.split("//")[-1].split("/", 1)
        if bucket.find("@") > 0:
            bucket, region = bucket.split("@")
        else:
            region = None
        return RemoteFile(bucket, key, region)
    else:
        raise NotImplementedError("Unexpected object store %s" % fname)

# ## Establish connection

def connect(fname):
    if fname.startswith("s3://"):
        return _connect_s3(fname)
    else:
        raise NotImplementedError("Unexpected object store %s" % fname)

def _connect_s3(fname=None):
    return boto.s3.connect_to_region(_get_region_s3(fname))

# ## Get region for object stores with multi-region support

def _get_region_s3(fname=None):
    """Retrieve region from standard environmental variables or file name.

    http://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-started.html#cli-environment
    """
    if fname:
        r_fname = parse_remote(fname)
        if r_fname.region:
            return r_fname.region
    return os.environ.get("AWS_DEFAULT_REGION", "us-east-1")

def default_region(fname):
    if fname.startswith("s3://"):
        return _get_region_s3()
    else:
        raise NotImplementedError("Unexpected object store %s" % fname)

# ## Retrieve files

def _s3_download_cl(fname):
    """Provide potentially streaming download from S3 using gof3r

    Selects the correct endpoint for non us-east support:
    http://docs.aws.amazon.com/general/latest/gr/rande.html#s3_region
    """
    r_fname = parse_remote(fname)
    cmd = ["gof3r", "get", "--no-md5", "-k", r_fname.key, "-b", r_fname.bucket]
    region = _get_region_s3(fname)
    if region != "us-east-1":
        cmd.append("--endpoint=s3-%s.amazonaws.com" % region)
    return cmd

def download(fname, input_dir, dl_dir=None):
    if fname.startswith("s3://"):
        r_fname = parse_remote(fname)
        if not dl_dir:
            dl_dir = utils.safe_makedir(os.path.join(input_dir, r_fname.bucket, os.path.dirname(r_fname.key)))
        out_file = os.path.join(dl_dir, os.path.basename(r_fname.key))
        if not utils.file_exists(out_file):
            with file_transaction({}, out_file) as tx_out_file:
                cmd = _s3_download_cl(fname) + ["-p", tx_out_file]
                subprocess.check_call(cmd)
        return out_file
    elif not fname or os.path.exists(fname):
        return fname
    else:
        raise NotImplementedError("Unexpected object store %s" % fname)

def cl_input(fname, unpack=True, anonpipe=True):
    """Return command line input for a file, handling streaming remote cases.
    """
    if not fname:
        return fname
    elif fname.startswith("s3://"):
        cmd = _s3_download_cl(fname)
        if fname.endswith(".gz") and unpack:
            cmd += " | gunzip -c"
        if anonpipe:
            cmd = "<(" + cmd + ")"
        return cmd
    else:
        return fname

# ## List buckets

def list(remote_dirname):
    if remote_dirname.startswith("s3://"):
        return _s3_list(remote_dirname)
    else:
        raise NotImplementedError("Unexpected object store %s" % remote_dirname)

def _s3_list(remote_dirname):
    r_fname = parse_remote(remote_dirname)
    conn = boto.connect_s3(remote_dirname)
    bucket = conn.get_bucket(r_fname.bucket)
    out = []
    region = "@%s" % r_fname.region if r_fname.region else ""
    for key in bucket.get_all_keys(prefix=r_fname.key):
        out.append("s3://%s%s/%s" % (r_fname.bucket, region, key.name))
    return out

# ## File handles

def open(fname):
    """Provide a handle-like object for streaming
    """
    if fname.startswith("s3://"):
        return _s3_handle(fname)
    else:
        raise NotImplementedError("Unexpected object store %s" % fname)

def _s3_handle(fname):
    """Return a handle like object for streaming from S3.
    """
    class S3Handle:
        def __init__(self, key):
            self._key = key
            self._iter = self._line_iter()
        def _line_iter(self):
            """From mrjob: https://github.com/Yelp/mrjob/blob/master/mrjob/util.py
            """
            buf = ""
            search_offset = 0
            for chunk in self._chunk_iter():
                buf += chunk
                start = 0
                while True:
                    end = buf.find("\n", start + search_offset) + 1
                    if end:  # if find() returned -1, end would be 0
                        yield buf[start:end]
                        start = end
                        # reset the search offset
                        search_offset = 0
                    else:
                        # this will happen eventually
                        buf = buf[start:]
                        # set search offset so we do not need to scan this part of the buffer again
                        search_offset = len(buf)
                        break
                if buf:
                    yield buf + '\n'
        def _chunk_iter(self):
            dec = zlib.decompressobj(16 | zlib.MAX_WBITS) if self._key.name.endswith(".gz") else None
            for chunk in self._key:
                if dec:
                    chunk = dec.decompress(chunk)
                if chunk:
                    yield chunk
        def __enter__(self):
            return self
        def __exit__(self, *args):
            self.close()
        def __iter__(self):
            return self
        def read(self, size):
            return self._key.read(size)
        def next(self):
            return self._iter.next()
        def close(self):
            self._key.close(fast=True)

    r_fname = parse_remote(fname)
    s3 = connect(fname)
    try:
        s3b = s3.get_bucket(r_fname.bucket)
    except boto.exception.S3ResponseError, e:
        # if we don't have bucket permissions but folder permissions, try without validation
        if e.status == 403:
            s3b = s3.get_bucket(r_fname.bucket, validate=False)
        else:
            raise
    s3key = s3b.get_key(r_fname.key)
    return S3Handle(s3key)
