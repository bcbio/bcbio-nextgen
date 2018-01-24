"""Extract files from processing run into output directory, organized by sample.
"""
import os
import shutil

from bcbio import utils
from bcbio.log import logger
from bcbio.upload import shared

def update_file(finfo, sample_info, config, pass_uptodate=False):
    """Update the file in local filesystem storage.
    """
    storage_dir = utils.safe_makedir(_get_storage_dir(finfo, config))
    if finfo.get("type") == "directory":
        return _copy_finfo_directory(finfo, storage_dir)
    else:
        return _copy_finfo(finfo, storage_dir, pass_uptodate=pass_uptodate)

def get_upload_path(finfo, sample_info, config):
    """"Dry" update the file: only return the upload path
    """
    try:
        storage_dir = _get_storage_dir(finfo, config)
    except ValueError:
        return None
    
    if finfo.get("type") == "directory":
        return _get_dir_upload_path(finfo, storage_dir)
    else:
        return _get_file_upload_path(finfo, storage_dir)

def _get_storage_dir(finfo, config):
    # skip if we have no directory to upload to
    if "dir" not in config:
        raise ValueError("Expect `dir` in upload specification: "
                         "http://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html#upload")
    if "run" in finfo:
        storage_dir = os.path.join(config["dir"], finfo["run"])
    elif "sample" in finfo:
        storage_dir = os.path.join(config["dir"], finfo["sample"])
    else:
        raise ValueError("Unexpected input file information: %s" % finfo)
    if "dir" in finfo:
        storage_dir = os.path.join(storage_dir, finfo["dir"])
    return storage_dir

def _get_file_upload_path(finfo, storage_dir):
    if "sample" in finfo and "ext" in finfo and "type" in finfo:
        out_file = os.path.join(storage_dir, "%s-%s%s%s" % (finfo["sample"], finfo["ext"],
                                                            "-" if (".txt" in finfo["type"]) else ".",
                                                            finfo["type"]))
    elif "batch" in finfo and "ext" in finfo and "type" in finfo:
        out_file = os.path.join(storage_dir, "%s-%s%s%s" % (finfo["batch"], finfo["ext"],
                                                            "-" if (".txt" in finfo["type"]) else ".",
                                                            finfo["type"]))
    else:
        out_file = os.path.join(storage_dir, os.path.basename(finfo["path"]))
    return os.path.abspath(out_file)

def _get_dir_upload_path(finfo, storage_dir):
    return os.path.abspath(os.path.join(storage_dir, finfo["ext"]))

def _copy_finfo(finfo, storage_dir, pass_uptodate=False):
    """Copy a file into the output storage directory.
    """
    out_file = _get_file_upload_path(finfo, storage_dir)
    if not shared.up_to_date(out_file, finfo):
        logger.info("Storing in local filesystem: %s" % out_file)
        shutil.copy(finfo["path"], out_file)
        return out_file
    if pass_uptodate:
        return out_file

def _copy_finfo_directory(finfo, out_dir):
    """Copy a directory into the final output directory.
    """
    out_dir = _get_dir_upload_path(finfo, out_dir)
    if not shared.up_to_date(out_dir, finfo):
        logger.info("Storing directory in local filesystem: %s" % out_dir)
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        shutil.copytree(finfo["path"], out_dir)
        for tmpdir in ["tx", "tmp"]:
            if os.path.exists(os.path.join(out_dir, tmpdir)):
                shutil.rmtree(os.path.join(out_dir, tmpdir))
        os.utime(out_dir, None)
    return out_dir
