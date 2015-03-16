"""Manage pushing and pulling files from an object store like Amazon Web Services S3.
"""
import os

import boto

from bcbio import utils

# ## definitions

SUPPORTED_REMOTES = ("s3://",)
BIODATA_INFO = \
  {"S3": {"bucket": "biodata",
          "key": "prepped/{build}/{build}-{target}.tar.gz"}}

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

# ## Establish connection

def connect(fname):
    if fname.startswith("s3://"):
        return _connect_s3()
    else:
        raise NotImplementedError("Unexpected object store %s" % fname)

def _connect_s3():
    return boto.s3.connect_to_region(_get_region_s3())

# ## Get region for object stores with multi-region support

def _get_region_s3():
    """Retrieve region from standard environmental variables.

    http://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-started.html#cli-environment
    """
    return os.environ.get("AWS_DEFAULT_REGION", "us-east-1")

# ## Retrieve file

# ## File handle
