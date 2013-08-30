"""Extract files from processing run into output directory, organized by sample.
"""
import os
import shutil

from bcbio import utils
from bcbio.log import logger
from bcbio.upload import shared

def copy_finfo(finfo, storage_dir):
    if finfo.has_key("sample"):
        out_file = os.path.join(storage_dir, "%s-%s.%s" % (finfo["sample"], finfo["ext"],
                                                           finfo["type"]))
    else:
        out_file = os.path.join(storage_dir, os.path.basename(finfo["path"]))
    out_file = os.path.abspath(out_file)
    if not shared.up_to_date(out_file, finfo):
        logger.info("Storing in local filesystem: %s" % out_file)
        shutil.copy(finfo["path"], out_file)
        return out_file

def update_file(finfo, sample_info, config):
    """Update the file in local filesystem storage.
    """
    # skip if we have no directory to upload to
    if "dir" not in config:
        return
    if finfo.has_key("sample"):
        storage_dir = utils.safe_makedir(os.path.join(config["dir"], finfo["sample"]))
    elif finfo.has_key("run"):
        storage_dir = utils.safe_makedir(os.path.join(config["dir"], finfo["run"]))
    else:
        raise ValueError("Unexpected input file information: %s" % finfo)
    return copy_finfo(finfo, storage_dir)
