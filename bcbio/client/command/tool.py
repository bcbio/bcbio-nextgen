"""Tools and utilities commands."""
import abc
import os
import sys
import contextlib

try:
    import progressbar
except ImportError:
    pass
import requests
from bcbio.client import base
from bcbio.distributed import objectstore


class Downloader(base.Command):

    """Base class for all the downloaders."""

    def __init__(self, parent, parser):
        super(Downloader, self).__init__(parent, parser)
        self._widgets = None
        self._progress = 0
        self._progress_bar = None
        self._size = 0

    def _update(self, chunk_size):
        """Update the progress bar."""
        if self.args.file != "-":
            self._progress += chunk_size
            self._progress_bar.update(self._progress)

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

    def prologue(self):
        """Executed once before the command running."""
        if self.args.file != "-":
            self._widgets = [
                os.path.basename(self.args.file), ": ",
                progressbar.Bar(marker="#", left="[", right="|"),
                progressbar.Percentage(), " ]",
                progressbar.FileTransferSpeed(), " ",
                progressbar.ETA()
            ]
            self._progress_bar = progressbar.ProgressBar(widgets=self._widgets,
                                                         maxval=self._size)
            self._progress_bar.start()

    def epilogue(self):
        """Executed once after the command running."""
        if self.args.file != "-":
            self._progress_bar.finish()


class S3Downloader(Downloader):

    """Download file from Amazon Simple Storage Service."""

    def __init__(self, parent, parser):
        super(S3Downloader, self).__init__(parent, parser)
        self._s3_file = None
        self._handle = None

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
        region = "-%s" % self.args.region if self.args.region else ""

        if "AWS_ACCESS_KEY_ID" not in os.environ:
            s3_file = "https://s3{region}.amazonaws.com/{bucket}/{key}"
            s3_file = s3_file.format(region=region, key=self.args.key,
                                     bucket=self.args.bucket)
            self._handle = requests.get(s3_file, stream=True)
            self._size = int(self._handle.headers['Content-Length'].strip())

        else:
            s3_file = "s3://{bucket}{region}{key}"
            s3_file = s3_file.format(region=self.args.region,
                                     key=self.args.key,
                                     bucket=self.args.bucket)
            self._handle = objectstore.AmazonS3.open(s3_file)
            self._size = self._handle.size

        super(S3Downloader, self).prologue()

    def work(self):
        """Run the command with the received information."""
        if "AWS_ACCESS_KEY_ID" not in os.environ:
            iterator = self._handle.iter_content(chunk_size=10240)
        else:
            iterator = self._handle

        with self._get_handler() as output:
            for chunk in iterator:
                output.write(chunk)
                self._update(len(chunk))


class BlobDownloader(Downloader):

    """Download file from Azure Blob storage service."""

    def __init__(self, parent, parser):
        super(BlobDownloader, self).__init__(parent, parser)
        self._service = None
        self._handle = None

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
        service = objectstore.AzureBlob.connect(file_info)
        self._handle = objectstore.BlobHandle(blob_service=service,
                                              container=self.args.container,
                                              blob=self.args.blob,
                                              chunk_size=10240)
        self._size = self._handle.blob_properties.get('content-length')
        super(BlobDownloader, self).prologue()

    def work(self):
        """Run the command with the received information."""
        with self._get_handler() as output:
            for chunk in self._handle:
                output.write(chunk)
                self._update(10240)
