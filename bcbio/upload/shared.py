"""Shared functionality for managing upload of final files.
"""
import os
import datetime

from bcbio import utils

def get_file_timestamp(f):
    return datetime.datetime.fromtimestamp(os.path.getmtime(f))

def up_to_date(new, orig):
    return (utils.file_exists(new) and
            get_file_timestamp(new) >= orig["mtime"])
