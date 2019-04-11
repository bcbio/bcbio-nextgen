"""
Manage pushing and pulling files from an object store like
Amazon Web Services S3.
"""
# pylint: disable=redefined-builtin

import abc
import collections
import os
import re
import subprocess
import sys
import time
import zlib

import six

from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio import utils

SUPPORTED_REMOTES = ("s3://",)
BIODATA_INFO = {"s3": "s3://biodata/prepped/{build}/{build}-{target}.tar.gz"}
REGIONS_NEWPERMS = {"s3": ["eu-central-1"]}

@six.add_metaclass(abc.ABCMeta)
class FileHandle(object):

    """Contract class for the file handle."""

    def __init__(self):
        self._iter = self._line_iter()

    def __enter__(self):
        """Define what the context manager should do at the beginning
        of the block created by the with statement.
        """
        return self

    def __exit__(self, *args):
        """Define what the context manager should do after its block
        has been executed (or terminates).
        """
        self.close()

    def __iter__(self):
        """Return the iterator for the current file."""
        return self

    def _line_iter(self):
        """Storage manager file iterator splits by buffer size instead
        of by newline. This wrapper puts them back into lines.

        From mrjob: https://github.com/Yelp/mrjob/blob/master/mrjob/util.py
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
                    # set search offset so we do not need to scan this part
                    # of the buffer again
                    search_offset = len(buf)
                    break
        if buf:
            yield buf + '\n'

    @abc.abstractmethod
    def _chunk_iter(self):
        """Chunk iterator over the received file."""
        pass

    @abc.abstractmethod
    def read(self, size=sys.maxsize):
        """Read at most size bytes from the file (less if the read hits EOF
        before obtaining size bytes).
        """
        pass

    @abc.abstractmethod
    def next(self):
        """Return the next item from the container."""
        pass

    @abc.abstractmethod
    def close(self):
        """Close the file handle."""
        pass


class S3Handle(FileHandle):

    """File object for the Amazon S3 files."""

    def __init__(self, key):
        super(S3Handle, self).__init__()
        self._key = key
        if self._key.name.endswith(".gz"):
            decompress = zlib.decompressobj(16 | zlib.MAX_WBITS)
            self._decompress = decompress.decompress
        else:
            self._decompress = lambda value: value

    def _chunk_iter(self):
        """Iterator over the S3 file."""
        for chunk in self._key:
            yield self._decompress(chunk)

    def read(self, size=sys.maxsize):
        """Read at most size bytes from the file (less if the read hits EOF
        before obtaining size bytes).
        """
        return self._key.read(size)

    def next(self):
        """Return the next item from the container."""
        return next(self._iter)

    def close(self):
        """Close the file handle."""
        self._key.close(fast=True)


class BlobHandle(FileHandle):

    """File object for the Azure Blob files."""

    def __init__(self, blob_service, container, blob, chunk_size):
        super(BlobHandle, self).__init__()
        self._blob_service = blob_service
        self._container_name = container
        self._blob_name = blob
        self._chunk_size = chunk_size
        self._blob_properties = {}
        self._pointer = 0

        if blob.endswith(".gz"):
            decompress = zlib.decompressobj(16 | zlib.MAX_WBITS)
            self._decompress = decompress.decompress
        else:
            self._decompress = lambda value: value

    @property
    def blob_properties(self):
        """Returns all user-defined metadata, standard HTTP properties,
        and system properties for the blob.
        """
        if not self._blob_properties:
            self._blob_properties = self._blob_service.get_blob_properties(
                container_name=self._container_name,
                blob_name=self._blob_name)
        return self._blob_properties

    def _chunk_offsets(self):
        """Iterator over chunk offests."""
        index = 0
        blob_size = self.blob_properties.get('content-length')
        while index < blob_size:
            yield index
            index = index + self._chunk_size

    def _chunk_iter(self):
        """Iterator over the blob file."""
        for chunk_offset in self._chunk_offsets():
            yield self._download_chunk(chunk_offset=chunk_offset,
                                       chunk_size=self._chunk_size)

    def _download_chunk_with_retries(self, chunk_offset, chunk_size,
                                     retries=3, retry_wait=1):
        """Reads or downloads the received blob from the system."""
        import azure
        while True:
            try:
                chunk = self._download_chunk(chunk_offset, chunk_size)
            except azure.WindowsAzureError:
                if retries > 0:
                    retries = retries - 1
                    time.sleep(retry_wait)
                else:
                    raise
            else:
                return chunk

    def _download_chunk(self, chunk_offset, chunk_size):
        """Reads or downloads the received blob from the system."""
        range_id = 'bytes={0}-{1}'.format(
            chunk_offset, chunk_offset + chunk_size - 1)

        return self._blob_service.get_blob(
            container_name=self._container_name,
            blob_name=self._blob_name,
            x_ms_range=range_id)

    def read(self, size=sys.maxsize):
        """Read at most size bytes from the file (less if the read hits EOF
        before obtaining size bytes).
        """
        blob_size = int(self.blob_properties.get('content-length'))
        if self._pointer < blob_size:
            chunk = self._download_chunk_with_retries(
                chunk_offset=self._pointer, chunk_size=size)
            self._pointer += size
            return chunk

    def next(self):
        """Return the next item from the container."""
        return next(self._iter)

    def close(self):
        """Close the file handle."""
        pass


@six.add_metaclass(abc.ABCMeta)
class StorageManager(object):

    """The contract class for all the storage managers."""

    @abc.abstractmethod
    def check_resource(self, resource):
        """Check if the received resource can be processed by
        the current storage manager.
        """
        pass

    @abc.abstractmethod
    def parse_remote(self, filename):
        """Parse a remote filename in order to obtain information
        related to received resource.
        """
        pass

    @abc.abstractmethod
    def connect(self, resource):
        """Return a connection object pointing to the endpoint
        associated to the received resource.
        """
        pass

    @abc.abstractmethod
    def download(self, filename, input_dir, dl_dir=None):
        """Download the resource from the storage."""
        pass

    @abc.abstractmethod
    def list(self, path):
        """Return a list containing the names of the entries in the directory
        given by path. The list is in arbitrary order.
        """
        pass

    @abc.abstractmethod
    def open(self, filename):
        """Provide a handle-like object for streaming."""
        pass


class AmazonS3(StorageManager):

    """Amazon Simple Storage Service (Amazon S3) Manager."""

    _DEFAULT_REGION = "us-east-1"
    _REMOTE_FILE = collections.namedtuple(
        "RemoteFile", ["store", "bucket", "key", "region"])
    _S3_FILE = "s3://%(bucket)s%(region)s/%(key)s"

    def __init__(self):
        super(AmazonS3, self).__init__()

    @classmethod
    def parse_remote(cls, filename):
        """Parses a remote filename into bucket and key information.

        Handles S3 with optional region name specified in key:
            BUCKETNAME@REGIONNAME/KEY
        """
        parts = filename.split("//")[-1].split("/", 1)
        bucket, key = parts if len(parts) == 2 else (parts[0], None)
        if bucket.find("@") > 0:
            bucket, region = bucket.split("@")
        else:
            region = None

        return cls._REMOTE_FILE("s3", bucket, key, region)

    @classmethod
    def _cl_aws_cli(cls, file_info, region):
        """Command line required for download using the standard AWS
        command line interface.
        """
        s3file = cls._S3_FILE % {"bucket": file_info.bucket,
                                 "key": file_info.key,
                                 "region": ""}
        command = [os.path.join(os.path.dirname(sys.executable), "aws"),
                   "s3", "cp", "--region", region, s3file]
        return (command, "awscli")

    @staticmethod
    def _cl_gof3r(file_info, region):
        """Command line required for download using gof3r."""
        command = ["gof3r", "get", "--no-md5",
                   "-k", file_info.key,
                   "-b", file_info.bucket]
        if region != "us-east-1":
            command += ["--endpoint=s3-%s.amazonaws.com" % region]
        return (command, "gof3r")

    @classmethod
    def _download_cl(cls, filename):
        """Provide potentially streaming download from S3 using gof3r
        or the AWS CLI.

        Selects the correct endpoint for non us-east support:
            http://docs.aws.amazon.com/general/latest/gr/rande.html#s3_region

        In eu-central-1 gof3r does not support new AWS signatures,
        so we fall back to the standard AWS commandline interface:
            https://github.com/rlmcpherson/s3gof3r/issues/45
        """
        file_info = cls.parse_remote(filename)
        region = cls.get_region(filename)
        if region in REGIONS_NEWPERMS["s3"]:
            return cls._cl_aws_cli(file_info, region)
        else:
            return cls._cl_gof3r(file_info, region)

    @classmethod
    def get_region(cls, resource=None):
        """Retrieve region from standard environmental variables
        or file name.

        More information of the following link: http://goo.gl/Vb9Jky
        """
        if resource:
            resource_info = cls.parse_remote(resource)
            if resource_info.region:
                return resource_info.region

        return os.environ.get("AWS_DEFAULT_REGION", cls._DEFAULT_REGION)

    @classmethod
    def check_resource(cls, resource):
        """Check if the received resource can be processed by
        the current storage manager.
        """
        if resource and resource.startswith("s3://"):
            return True

        return False

    @classmethod
    def connect(cls, resource):
        """Connect to this Region's endpoint.

        Returns a connection object pointing to the endpoint associated
        to the received resource.
        """
        import boto
        return boto.s3.connect_to_region(cls.get_region(resource))

    @classmethod
    def download(cls, filename, input_dir, dl_dir=None):
        """Provide potentially streaming download from S3 using gof3r
        or the AWS CLI.
        """
        file_info = cls.parse_remote(filename)
        if not dl_dir:
            dl_dir = os.path.join(input_dir, file_info.bucket,
                                  os.path.dirname(file_info.key))
            utils.safe_makedir(dl_dir)

        out_file = os.path.join(dl_dir, os.path.basename(file_info.key))

        if not utils.file_exists(out_file):
            with file_transaction({}, out_file) as tx_out_file:
                command, prog = cls._download_cl(filename)
                if prog == "gof3r":
                    command.extend(["-p", tx_out_file])
                elif prog == "awscli":
                    command.extend([tx_out_file])
                else:
                    raise NotImplementedError(
                        "Unexpected download program %s" % prog)
                subprocess.check_call(command)
        return out_file

    @classmethod
    def cl_input(cls, filename, unpack=True, anonpipe=True):
        """Return command line input for a file, handling streaming
        remote cases.
        """
        command, prog = cls._download_cl(filename)
        if prog == "awscli":
            command.append("-")

        command = " ".join(command)
        if filename.endswith(".gz") and unpack:
            command = "%(command)s | gunzip -c" % {"command": command}
        elif filename.endswith(".bz2") and unpack:
            command = "%(command)s | bunzip2 -c" % {"command": command}
        if anonpipe:
            command = "<(%(command)s)" % {"command": command}

        return command

    @classmethod
    def list(cls, path):
        """Return a list containing the names of the entries in the directory
        given by path. The list is in arbitrary order.
        """
        file_info = cls.parse_remote(path)
        connection = cls.connect(path)
        bucket = connection.get_bucket(file_info.bucket)
        region = "@%s" % file_info.region if file_info.region else ""

        output = []
        for key in bucket.get_all_keys(prefix=file_info.key):
            output.append(cls._S3_FILE % {"bucket": file_info.bucket,
                                          "key": key.name,
                                          "region": region})
        return output

    @classmethod
    def open(cls, filename):
        """Return a handle like object for streaming from S3."""
        import boto
        file_info = cls.parse_remote(filename)
        connection = cls.connect(filename)
        try:
            s3_bucket = connection.get_bucket(file_info.bucket)
        except boto.exception.S3ResponseError as error:
            # if we don't have bucket permissions but folder permissions,
            # try without validation
            if error.status == 403:
                s3_bucket = connection.get_bucket(file_info.bucket,
                                                  validate=False)
            else:
                raise

        s3_key = s3_bucket.get_key(file_info.key)
        if s3_key is None:
            raise ValueError("Did not find S3 key: %s" % filename)
        return S3Handle(s3_key)


class AzureBlob(StorageManager):

    """Azure Blob storage service manager."""

    _BLOB_FILE = ("https://{storage}.blob.core.windows.net/"
                  "{container}/{blob}")
    _REMOTE_FILE = collections.namedtuple(
        "RemoteFile", ["store", "storage", "container", "blob"])
    _URL_FORMAT = re.compile(r'http.*\/\/(?P<storage>[^.]+).blob.core.windows.net\/'
                             r'(?P<container>[^/]+)\/*(?P<blob>.*)')
    _BLOB_CHUNK_DATA_SIZE = 4 * 1024 * 1024

    def __init__(self):
        super(AzureBlob, self).__init__()

    @classmethod
    def check_resource(cls, resource):
        """Check if the received resource can be processed by
        the current storage manager.
        """
        return cls._URL_FORMAT.match(resource or "")

    @classmethod
    def parse_remote(cls, filename):
        """Parses a remote filename into blob information."""
        blob_file = cls._URL_FORMAT.search(filename)
        return cls._REMOTE_FILE("blob",
                                storage=blob_file.group("storage"),
                                container=blob_file.group("container"),
                                blob=blob_file.group("blob"))

    @classmethod
    def connect(cls, resource):
        """Returns a connection object pointing to the endpoint
        associated to the received resource.
        """
        from azure import storage as azure_storage
        file_info = cls.parse_remote(resource)
        return azure_storage.BlobService(file_info.storage)

    @classmethod
    def download(cls, filename, input_dir, dl_dir=None):
        """Download the resource from the storage."""
        file_info = cls.parse_remote(filename)
        if not dl_dir:
            dl_dir = os.path.join(input_dir, file_info.container,
                                  os.path.dirname(file_info.blob))
            utils.safe_makedir(dl_dir)

        out_file = os.path.join(dl_dir, os.path.basename(file_info.blob))

        if not utils.file_exists(out_file):
            with file_transaction({}, out_file) as tx_out_file:
                blob_service = cls.connect(filename)
                blob_service.get_blob_to_path(
                    container_name=file_info.container,
                    blob_name=file_info.blob,
                    file_path=tx_out_file)
        return out_file

    @classmethod
    def list(cls, path):
        """Return a list containing the names of the entries in the directory
        given by path. The list is in arbitrary order.
        """
        import azure
        from azure import storage as azure_storage
        output = []
        path_info = cls.parse_remote(path)
        blob_service = azure_storage.BlobService(path_info.storage)
        try:
            blob_enum = blob_service.list_blobs(path_info.container)
        except azure.WindowsAzureMissingResourceError:
            return output

        for item in blob_enum:
            output.append(cls._BLOB_FILE.format(storage=path_info.storage,
                                                container=path_info.container,
                                                blob=item.name))
        return output

    @classmethod
    def open(cls, filename):
        """Provide a handle-like object for streaming."""
        file_info = cls.parse_remote(filename)
        blob_service = cls.connect(filename)
        return BlobHandle(blob_service=blob_service,
                          container=file_info.container,
                          blob=file_info.blob,
                          chunk_size=cls._BLOB_CHUNK_DATA_SIZE)

class ArvadosKeep:
    """Files stored in Arvados Keep. Partial implementation, integration in bcbio-vm.
    """
    @classmethod
    def check_resource(self, resource):
        return resource.startswith("keep:")
    @classmethod
    def download(self, filename, input_dir, dl_dir=None):
        return None

class SevenBridges:
    """Files stored in SevenBridges. Partial implementation, integration in bcbio-vm.
    """
    @classmethod
    def check_resource(self, resource):
        return resource.startswith("sbg:")
    @classmethod
    def download(self, filename, input_dir, dl_dir=None):
        return None

class DNAnexus:
    """Files stored in DNAnexus. Partial implementation, integration in bcbio-vm.
    """
    @classmethod
    def check_resource(self, resource):
        return resource.startswith("dx:")
    @classmethod
    def download(self, filename, input_dir, dl_dir=None):
        return None

class GoogleCloud:
    """Files stored in Google Cloud Storage. Partial implementation, integration in bcbio-vm.
    """
    @classmethod
    def check_resource(self, resource):
        return resource.startswith("gs:")
    @classmethod
    def download(self, filename, input_dir, dl_dir=None):
        return None

class RegularServer:
    """Files stored in FTP/http that can be downloaded by wget
    """

    @classmethod
    def _parse_url(self, fn):
        regex = re.compile(r'^(?:http|ftp)s?://.*(.fastq.gz)[^/]*|' # http:// or https://
                           r'^(?:http|ftp)s?://.*(.fastq)[^/]*|'
                           r'^(?:http|ftp)s?://.*(.bam)[^/]*'
                           , re.IGNORECASE)
        return regex.match(fn)

    @classmethod
    def check_resource(self, resource):
        if self._parse_url(resource):
            return True

    @classmethod
    def download(self, filename, input_dir, dl_dir=None):
        match = self._parse_url(filename)
        file_info = os.path.basename(filename)
        ext = file_info.find(match.group(1))
        name = file_info[:ext]
        if not dl_dir:
            dl_dir = os.path.join(input_dir, name)
            utils.safe_makedir(dl_dir)

        fixed_name = "%s%s" % (name, match.group(1))
        out_file = os.path.join(dl_dir, fixed_name)

        if not utils.file_exists(out_file):
            with file_transaction({}, out_file) as tx_out_file:
                cmd = "wget -O {tx_out_file} {filename}"
                do.run(cmd.format(**locals()), "Download %s" % out_file)

        return out_file

def _get_storage_manager(resource):
    """Return a storage manager which can process this resource."""
    for manager in (AmazonS3, ArvadosKeep, SevenBridges, DNAnexus, AzureBlob, GoogleCloud, RegularServer):
        if manager.check_resource(resource):
            return manager()

    raise ValueError("Unexpected object store  %(resource)s" %
                     {"resource": resource})


def is_remote(fname):
    """Check if the received file is recognised by one of
    the available storage managers.
    """
    try:
        _get_storage_manager(fname)
    except ValueError:
        return False

    return True


def file_exists_or_remote(fname):
    """Check if a file exists or is accessible remotely."""
    return is_remote(fname) or utils.file_exists(fname)


def default_region(fname):
    """Return the default region for the received resource.

    Note:
        This feature is available only for AmazonS3 storage manager.
    """
    manager = _get_storage_manager(fname)
    if hasattr(manager, "get_region"):
        return manager.get_region()

    raise NotImplementedError("Unexpected object store %s" % fname)


def connect(filename):
    """Returns a connection object pointing to the endpoint associated
    to the received resource.
    """
    manager = _get_storage_manager(filename)
    return manager.connect(filename)


def download(fname, input_dir, dl_dir=None):
    """Download the resource from the storage."""
    try:
        manager = _get_storage_manager(fname)
    except ValueError:
        return fname
    return manager.download(fname, input_dir, dl_dir)


def cl_input(fname, unpack=True, anonpipe=True):
    """Return command line input for a file, handling streaming
    remote cases.
    """
    try:
        manager = _get_storage_manager(fname)
    except ValueError:
        return fname

    return manager.cl_input(fname, unpack, anonpipe)


def list(remote_dirname):
    """Return a list containing the names of the entries in the directory
    given by path. The list is in arbitrary order.
    """
    manager = _get_storage_manager(remote_dirname)
    return manager.list(remote_dirname)


def open_file(fname):
    """Provide a handle-like object for streaming."""
    manager = _get_storage_manager(fname)
    return manager.open(fname)


def parse_remote(fname):
    """Parses a remote filename in order to obtain information
    related to received resource.
    """
    manager = _get_storage_manager(fname)
    return manager.parse_remote(fname)
