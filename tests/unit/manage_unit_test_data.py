"""
Manages the unit test data. Run with no arguments it checks to see if
the unit test data is stale and if it is updates it. Can also be used
to upload data to S3.

Note: for some reason the md5 hash for the tarred data directory is always
different, so there is no checking for staleness on upload.

"""

import subprocess
import os
import urllib2
import posixpath
import hashlib
import argparse
from boto.s3.connection import S3Connection
from boto.s3.key import Key
from boto.exception import S3CreateError


BUCKET = "bcbio_nextgen"
DATA_FILE = "unit_test_data.tar.gz"
DATA_URL = "https://s3.amazonaws.com/" + BUCKET + "/" + DATA_FILE
UNIT_TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def upload_unitdata_to_s3(access_key, secret_key, bucket_name, filename):
    print "Connecting to %s." % (bucket_name)
    conn = S3Connection(access_key, secret_key)
    # make the bucket if it doesn't already exist
    make_bucket(conn, bucket_name)
    bucket = conn.get_bucket(bucket_name)
    print "Uploading %s to %s under key %s." % (filename, bucket_name,
                                                DATA_FILE)
    upload_file_to_bucket(bucket, DATA_FILE, filename)
    # store a md5 sum as a file as well. might be possible to do
    # this with metadata too
    hash_string = md5sum(filename)
    print "Uploading %s as hash of %s under key %s." % (hash_string,
                                                        filename,
                                                        DATA_FILE + ".md5")
    upload_string_to_bucket(bucket, DATA_FILE + ".md5", md5sum(filename))
    print "Done!"


def upload_file_to_bucket(bucket, key_name, filename, public=True):
    k = Key(bucket)
    k.key = key_name
    k.set_contents_from_filename(filename)
    if public:
        k.set_acl('public-read')


def upload_string_to_bucket(bucket, key_name, s, public=True):
    k = Key(bucket)
    k.key = key_name
    k.set_contents_from_string(s)
    if public:
        k.set_acl('public-read')


def make_bucket(conn, bucket_name):
    try:
        bucket = conn.create_bucket(bucket_name)
    # expected if we are not the owner and someone else has made it
    except S3CreateError:
        pass


def tar_data_directory():
    print "Creating archive of unit test data directory."
    unit_test_file = posixpath.basename(DATA_URL)
    if os.path.exists(unit_test_file):
        os.remove(unit_test_file)
    cmd = ["tar", "-czvf", unit_test_file, UNIT_TEST_DATA_DIR]
    subprocess.check_call(cmd)
    return unit_test_file


def md5sum(filename):
    with open(filename, mode='rb') as f:
        d = hashlib.md5(f.read())
    return d.hexdigest()


def needs_update(filename, url):
    return not get_local_hash(filename) == get_remote_hash(url)


def get_local_hash(filename):
    return md5sum(filename)


def get_remote_hash(url):
    response = urllib2.urlopen(url)
    return response.read().strip()


def install_test_files():
    """Download required sequence and reference files.
    """
    local_data_file = posixpath.basename(DATA_URL)
    if os.path.exists(local_data_file):
        if not needs_update(local_data_file, DATA_URL + ".md5"):
            print "Unit test data already up to date."
            exit(1)
        else:
            print "Stale unit test data detected. Grabbing new data."
            os.unlink(local_data_file)

    response = urllib2.urlopen(DATA_URL)
    with open(local_data_file, "wb") as out_handle:
        out_handle.write(response.read())
    subprocess.check_call(["tar", "-zxvf", local_data_file])
    print "Done!"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--upload", help="upload data to S3",
                        action="store_true")
    parser.add_argument("--secret-key", help="secret key for uploading to S3.")
    parser.add_argument("--access-key", help="access key for uploading to S3.")
    args = parser.parse_args()
    if args.upload:
        if not (args.secret_key and args.access_key):
            parser.print_help()
            exit(1)
        else:
            data_file = tar_data_directory()
            upload_unitdata_to_s3(args.access_key, args.secret_key, BUCKET,
                                  data_file)

    else:
        install_test_files()
