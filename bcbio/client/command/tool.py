"""Tools and utilities commands."""
import abc
import os
import sys
import contextlib

import requests

from bcbio.client import base
from bcbio.distributed import objectstore


class Downloader(base.Command):

    """Base class for all the downloaders."""

    @contextlib.contextmanager
    def _get_handler(self):
        """Get the output file/stream handler."""
        if self.args.file and self.args.file != '-':
            file_handle = open(self.args.file, 'w')
        else:
            file_handle = sys.stdout

        try:
            yield file_handle
        finally:
            if file_handle is not sys.stdout:
                file_handle.close()

    @abc.abstractmethod
    def setup(self):
        """Extend the parser configuration in order to expose this command."""
        pass

    @abc.abstractmethod
    def work(self):
        """Run the command with the received information."""
        pass


class S3Downloader(Downloader):

    """Download file from Amazon Simple Storage Service."""

    def __init__(self, parent, parser):
        super(S3Downloader, self).__init__(parent, parser)
        self._s3_file = None

    def setup(self):
        """Extend the parser configuration in order to expose this command."""
        parser = self._parser.add_parser(
            "s3", help="Amazon Simple Storage Service (Amazon S3)")
        parser.add_argument(
            "--file", default="-",
            help="The path where the file should be copied.")
        parser.add_argument(
            "-k", "--key", required=True,
            help="The name of the file.")
        parser.add_argument(
            "-b", "--bucket", required=True,
            help="The name of the container that contains the file.")
        parser.add_argument(
            "-r", "--region", default="", help="The name of the region.")
        parser.set_defaults(work=self.run)

    def prologue(self):
        """Executed once before the command running."""
        s3_file = "https://s3{region}.amazonaws.com/{bucket}/{key}"
        region = "-%s" % self.args.region if self.args.region else ""
        self._s3_file = s3_file.format(region=region, bucket=self.args.bucket,
                                       key=self.args.key)

    def work(self):
        """Run the command with the received information."""
        response = requests.get(self._s3_file)
        with self._get_handler() as output:
            for chunk in response.iter_content(chunk_size=1024):
                output.write(chunk)


class BlobDownloader(Downloader):

    """Download file from Azure Blob storage service."""

    def __init__(self, parent, parser):
        super(BlobDownloader, self).__init__(parent, parser)
        self._service = None

    def setup(self):
        """Extend the parser configuration in order to expose this command."""
        parser = self._parser.add_parser(
            "blob", help="Azure Blob storage service")
        parser.add_argument(
            "--file", default="-",
            help="The path where the file should be copied.")
        parser.add_argument(
            "-b", "--blob", required=True,
            help="The name of the blob.")
        parser.add_argument(
            "-c", "--container", required=True,
            help="The name of the container that contains the blob. All "
                 "blobs must be in a container.")
        parser.add_argument(
            "-a", "--account_name",
            default=os.environ.get("STORAGE_ACCOUNT", None),
            help="The storage account name. All access to Azure Storage"
                 " is done through a storage account.")

    def prologue(self):
        """Executed once before the command running."""
        file_info = objectstore.AzureBlob.REMOTE_FILE(
            "blob", storage=self.args.account_name, blob=self.args.blob,
            container=self.args.container)
        self._service = objectstore.AzureBlob.connect(file_info)

    def work(self):
        """Run the command with the received information."""
        file_handle = objectstore.BlobHandle(blob_service=self._service,
                                             container=self.args.container,
                                             blob=self.args.blob,
                                             chunk_size=1024)

        with self._get_handler() as output:
            for chunk in file_handle:
                output.write(chunk)
