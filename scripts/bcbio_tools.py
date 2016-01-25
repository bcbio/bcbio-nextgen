#!/usr/bin/env python -E
"""bcbio-tools command line application."""

from __future__ import print_function
import argparse
import sys

from bcbio.client import base
from bcbio.client import groups


class BcbioTools(base.Client):

    """bcbio-tools command line application."""

    commands = [
        (groups.Downloader, "tools")
    ]

    def setup(self):
        """Extend the parser configuration in order to expose all
        the received commands.
        """
        self._parser = argparse.ArgumentParser(
            description=("Tools and utilities."))
        tools = self._parser.add_subparsers(title="[tools]")
        self._register_parser("tools", tools)


def main():
    """Run the bcbio-nextgen-vm command line application."""
    bcbio = BcbioTools(sys.argv[1:])
    bcbio.run()


if __name__ == "__main__":
    main()
