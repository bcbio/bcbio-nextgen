"""The commands used by the command line parser."""

from bcbio.client import base
from bcbio.client import command


class Downloader(base.Group):

    """Download file from a storage service."""

    commands = [
        (command.tool.S3Downloader, "storage_service"),
        (command.tool.BlobDownloader, "storage_service"),
    ]

    def setup(self):
        """Extend the parser configuration in order to expose this command."""
        download = self._parser.add_parser(
            "download", help="Download file to a storage service.")
        storage_service = download.add_subparsers(title="[storage manager]")

        self._register_parser("storage_service", storage_service)
