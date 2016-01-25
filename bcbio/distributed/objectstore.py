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

try:
    import azure
    from azure import storage as azure_storage
except ImportError:
    azure, azure_storage = None, None
import boto
import six
import requests

from bcbio.distributed.transaction import file_transaction
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
    def read(self, size):
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

    @property
    def size(self):
        """Return the file size."""
        return self._key.size

    def _chunk_iter(self):
        """Iterator over the S3 file."""
        for chunk in self._key:
            yield self._decompress(chunk)

    def read(self, size):
        """Read at most size bytes from the file (less if the read hits EOF
        before obtaining size bytes).
        """
        return self._key.read(size)

    def next(self):
        """Return the next item from the container."""
        return self._iter.next()

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

    def read(self, size):
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
        return self._iter.next()

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
    def connect(self, resource, context):
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

    @abc.abstractmethod
    def cl_input(self, filename, unpack=True, anonpipe=True):
        """Return command line input for a file, handling streaming
        remote cases.
        """
        pass

    @abc.abstractmethod
    def resource_exists(self, resource, context=None):
        """Check if the received resource exists on storage service."""
        pass

    @abc.abstractmethod
    def exists(self, container, filename, context=None):
        """Check if the received file name exists in the container."""
        pass

    @abc.abstractmethod
    def upload(self, path, filename, container, context=None):
        """Upload the received file.

        :path:      The path of the file that should be uploaded.
        :container: The name of the container.
        :filename:  The name of the item from the container.
        :context:   More information required by the storage manager.
        """
        pass


class AmazonS3(StorageManager):

    """Amazon Simple Storage Service (Amazon S3) Manager."""

    REMOTE_FILE = collections.namedtuple(
        "RemoteFile", ["store", "bucket", "key", "region"])

    _DEFAULT_REGION = "us-east-1"
    _S3_FILE = "s3://%(bucket)s%(region)s/%(key)s"
    _S3_URL = "https://s3{region}.amazonaws.com/{bucket}/{key}"
    _UPLOAD_HEADERS = {
        "x-amz-storage-class": "REDUCED_REDUNDANCY",
        "x-amz-server-side-encryption": "AES256",
    }

    def __init__(self):
        super(AmazonS3, self).__init__()

    @classmethod
    def get_bucket(cls, bucket_name):
        """Retrieves a bucket by name."""
        connection = boto.connect_s3()
        try:
            # If the bucket does not exist, an S3ResponseError
            # will be raised.
            bucket = connection.get_bucket(bucket_name)
        except boto.exception.S3ResponseError as exc:
            if exc.status == 404:
                bucket = connection.create_bucket(bucket_name)
            else:
                raise

        return bucket

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

        return cls.REMOTE_FILE("s3", bucket, key, region)

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

    @staticmethod
    def _cl_bcbio_tool(file_info, region):
        """Command line required for download file ussing bcbio_tools.py."""

        def _check_script(script):
            """Check if the script is available."""
            return subprocess.call(
                ["type", script], shell=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE) == 0

        binaries = ("bcbio_vm.py", "bcbio_tools.py")
        for script in binaries:
            if _check_script(script):
                command = [script, "download",
                           "--key", file_info.key,
                           "--bucket", file_info.bucket]
                return (command, script)

        raise ValueError("No commandline tool available.")

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

        if "AWS_ACCESS_KEY_ID" not in os.environ:
            return cls._cl_bcbio_tool(file_info, region)

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
    def connect(cls, resource, context=None):
        """Connect to this Region's endpoint.

        Returns a connection object pointing to the endpoint associated
        to the received resource.
        """
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
                elif prog == "bcbio_vm.py":
                    command[1:2] = ["tools", "download", "s3"]
                    command.extend(["--file", tx_out_file])
                elif prog == "bcbio_tools.py":
                    command.extend(["--file", tx_out_file])
                else:
                    raise NotImplementedError(
                        "Unexpected download program %s" % prog)
                subprocess.check_call(command, stdout=sys.stdout)
        return out_file

    @classmethod
    def cl_input(cls, filename, unpack=True, anonpipe=True):
        """Return command line input for a file, handling streaming
        remote cases.
        """
        command, prog = cls._download_cl(filename)
        if prog == "awscli":
            command.append("-")
        elif prog == "bcbio_vm.py":
            command[1:2] = ["tools", "download", "s3"]

        command = " ".join(command)
        if filename.endswith(".gz") and unpack:
            command = "%(command)s | gunzip -c" % {"command": command}
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
        file_info = cls.parse_remote(filename)
        s3_bucket = cls.get_bucket(file_info.bucket)
        s3_key = s3_bucket.get_key(file_info.key)
        if not s3_key:
            raise ValueError("The file %r is not available." %
                             file_info.key)
        return S3Handle(s3_key)

    @classmethod
    def resource_exists(cls, resource, context=None):
        """Check if the received key name exists in the bucket."""
        file_info = cls.parse_remote(resource)
        context = (context or {})["region"] = file_info.region

        return cls.exists(file_info.bucket, file_info.key, context)

    @classmethod
    def exists(cls, container, filename, context=None):
        """Check if the received key name exists in the bucket.

        :container: The name of the bucket.
        :filename:  The name of the key.
        :context:   More information required by the storage manager.
        """
        if "AWS_ACCESS_KEY_ID" not in os.environ:
            region = (context or {}).get("region", "")
            s3_url = cls._S3_URL.format(region=region, bucket=container,
                                        key=filename)
            response = requests.head(s3_url)
            # pylint: disable=no-member
            return response.status_code == requests.codes.ok

        bucket = cls.get_bucket(container)
        key = bucket.get_key(filename)
        return True if key else False

    @classmethod
    def upload(cls, path, filename, container, context=None):
        """Upload the received file.

        :path:      The path of the file that should be uploaded.
        :container: The name of the bucket.
        :filename:  The name of the key.
        :context:   More information required by the storage manager.
        """
        headers = (context or {}).get("headers", cls._UPLOAD_HEADERS)
        arguments = (context or {}).get("arguments", [])

        command = ["gof3r", "put", "-p", path,
                   "-k", filename, "-b", container]
        command.extend(arguments)

        if headers:
            for header, value in headers.items():
                command.extend(("-m", "{0}:{1}".format(header, value)))

        subprocess.check_call(command)


class AzureBlob(StorageManager):

    """Azure Blob storage service manager."""

    _BLOB_FILE = ("https://{storage}.blob.core.windows.net/"
                  "{container}/{blob}")
    REMOTE_FILE = collections.namedtuple(
        "RemoteFile", ["store", "storage", "container", "blob"])
    _URL_FORMAT = re.compile(r'http.*\/\/(?P<storage>[^.]+)[^/]+\/'
                             r'(?P<container>[^/]+)\/*(?P<blob>.*)')
    _BLOB_CHUNK_DATA_SIZE = 4 * 1024 * 1024

    def __init__(self):
        super(AzureBlob, self).__init__()

    @classmethod
    def _get_credentials(cls, context=None):
        """Get AzureBlob credentials from environment or context data."""
        credentials = (context or {}).get("credentials", {})

        if "storage_account" in credentials:
            account_name = credentials["storage_account"]
            account_key = credentials.get("storage_access_key", None)
        else:
            account_name = os.environ.get("STORAGE_ACCOUNT", None)
            account_key = os.environ.get("STORAGE_ACCESS_KEY", None)

        return (account_name, account_key)

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
        return cls.REMOTE_FILE("blob",
                               storage=blob_file.group("storage"),
                               container=blob_file.group("container"),
                               blob=blob_file.group("blob"))

    @classmethod
    def connect(cls, resource, context=None):
        """Returns a connection object pointing to the endpoint
        associated to the received resource.
        """
        if isinstance(resource, cls.REMOTE_FILE):
            file_info = resource
        else:
            file_info = cls.parse_remote(resource)

        account_name, account_key = cls._get_credentials(context)
        if account_name != file_info.storage and file_info.storage:
            account_name = file_info.storage
            account_key = None

        return azure_storage.BlobService(account_name=account_name,
                                         account_key=account_key)

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

    @classmethod
    def cl_input(cls, filename, unpack=True, anonpipe=True):
        """Return command line input for a file, handling streaming
        remote cases.
        """
        file_info = cls.parse_remote(filename)
        command = " ".join([
            "bcbio_tools.py", "download", "blob",
            "--blob", file_info.blob,
            "--container", file_info.container,
            "--account_name", file_info.storage,
        ])

        if filename.endswith(".gz") and unpack:
            command = "%(command)s | gunzip -c" % {"command": command}
        if anonpipe:
            command = "<(%(command)s)" % {"command": command}

        return command

    @classmethod
    def resource_exists(cls, resource, context=None):
        """Check if the received key name exists in the bucket."""
        file_info = cls.parse_remote(resource)
        return cls.exists(file_info.container, file_info.blob, context)

    @classmethod
    def exists(cls, container, filename, context=None):
        """Check if the received key name exists in the bucket.

        :container: The name of the container that contains the blob. All
                    blobs must be in a container.
        :filename:  The name of the blob.
        :context:   More information required by the storage manager.

        :notes:
            The context should contain the storage account name.
            All access to Azure Storage is done through a storage account.
        """

        account_name, _ = cls._get_credentials(context)
        if not account_name:
            raise ValueError("The account_name was not found in %s." %
                             {"container": "context: {0}".format(context)})

        blob_service = cls.connect(cls.REMOTE_FILE("blob", account_name,
                                                   container, filename))
        blob_handle = BlobHandle(blob_service=blob_service, blob=filename,
                                 container=container, chunk_size=32)
        try:
            blob_handle.read(1024)
        except azure.WindowsAzureMissingResourceError:
            return False

        return True

    @classmethod
    def upload(cls, path, filename, container, context=None):
        """Upload the received file.

        :path:       The path of the file that should be uploaded.
        :container:  The name of the container that contains the blob. All
                     blobs must be in a container.
        :filename:   The name of the blob.
        :context:    More information required by the storage manager.

        :notes:
            The context should contain the storage account name.
            All access to Azure Storage is done through a storage account.
        """
        file_info = cls.REMOTE_FILE("blob", "", container, filename)
        blob_service = cls.connect(file_info, context=context)
        # Ensure that the container exists.
        blob_service.create_container(container_name=container,
                                      x_ms_blob_public_access='container',
                                      fail_on_exist=False)
        # Upload the received file.
        blob_service.put_block_blob_from_path(container_name=container,
                                              blob_name=filename,
                                              file_path=path)


def _get_storage_manager(resource):
    """Return a storage manager which can process this resource."""
    for manager in (AmazonS3, AzureBlob):
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
    if is_remote(fname):
        return True
    else:
        return utils.file_exists(fname)


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


def open(fname):
    """Provide a handle-like object for streaming."""
    manager = _get_storage_manager(fname)
    return manager.open(fname)


def parse_remote(fname):
    """Parses a remote filename in order to obtain information
    related to received resource.
    """
    manager = _get_storage_manager(fname)
    return manager.parse_remote(fname)


def resource_exists(resource, context=None):
    """Check if the received resource exists on storage service."""
    manager = _get_storage_manager(resource)
    return manager.resource_exists(resource, context)
