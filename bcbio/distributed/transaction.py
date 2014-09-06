"""Handle file based transactions allowing safe restarts at any point.

To handle interrupts,this defines output files written to temporary
locations during processing and copied to the final location when finished.
This ensures output files will be complete independent of method of
interruption.
"""
import contextlib
import os
import uuid
import shutil
import tempfile

import toolz as tz

from bcbio import utils

@contextlib.contextmanager
def tx_tmpdir(data=None, base_dir=None, remove=True):
    """Context manager to create and remove a transactional temporary directory.

    Handles creating a transactional directory for running commands in. Will
    use either the current directory or a configured temporary directory.

    Creates an intermediary location and time specific directory for global
    temporary directories to prevent collisions.

    data can be the full world information object being process or a
    configuration dictionary.
    """
    if data and "config" in data:
        config_tmpdir = tz.get_in(("config", "resources", "tmp", "dir"), data)
    elif data:
        config_tmpdir = tz.get_in(("resources", "tmp", "dir"), data)
    else:
        config_tmpdir = None
    if config_tmpdir:
        tmp_dir_base = os.path.join(config_tmpdir, "bcbiotx", str(uuid.uuid1()))
    elif base_dir is not None:
        tmp_dir_base = os.path.join(base_dir, "tx")
    else:
        tmp_dir_base = os.path.join(os.getcwd(), "tx")
    utils.safe_makedir(tmp_dir_base)
    tmp_dir = tempfile.mkdtemp(dir=tmp_dir_base)
    utils.safe_makedir(tmp_dir)
    try:
        yield tmp_dir
    finally:
        if remove:
            try:
                shutil.rmtree(tmp_dir)
                if config_tmpdir:
                    shutil.rmtree(tmp_dir_base)
            except:
                pass

@contextlib.contextmanager
def file_transaction(*rollback_files):
    """Wrap file generation in a transaction, moving to output if finishes.
    """
    exts = {".vcf": ".idx", ".bam": ".bai", ".vcf.gz": ".tbi", ".bed.gz": ".tbi"}
    safe_names, orig_names = _flatten_plus_safe(rollback_files)
    _remove_files(safe_names)  # remove any half-finished transactions
    try:
        if len(safe_names) == 1:
            yield safe_names[0]
        else:
            yield tuple(safe_names)
    except:  # failure -- delete any temporary files
        _remove_files(safe_names)
        _remove_tmpdirs(safe_names)
        raise
    else:  # worked -- move the temporary files to permanent location
        for safe, orig in zip(safe_names, orig_names):
            if os.path.exists(safe):
                shutil.move(safe, orig)
                for check_ext, check_idx in exts.iteritems():
                    if safe.endswith(check_ext):
                        safe_idx = safe + check_idx
                        if os.path.exists(safe_idx):
                            shutil.move(safe_idx, orig + check_idx)
        _remove_tmpdirs(safe_names)

def _remove_tmpdirs(fnames):
    for x in fnames:
        xdir = os.path.dirname(os.path.abspath(x))
        if xdir and os.path.exists(xdir):
            shutil.rmtree(xdir, ignore_errors=True)

def _remove_files(fnames):
    for x in fnames:
        if x and os.path.exists(x):
            if os.path.isfile(x):
                os.remove(x)
            elif os.path.isdir(x):
                shutil.rmtree(x, ignore_errors=True)

def _flatten_plus_safe(rollback_files):
    """Flatten names of files and create temporary file names.
    """
    tx_files, orig_files = [], []
    for fnames in rollback_files:
        if isinstance(fnames, basestring):
            fnames = [fnames]
        for fname in fnames:
            basedir = utils.safe_makedir(os.path.join(os.path.dirname(fname), "tx"))
            tmpdir = utils.safe_makedir(tempfile.mkdtemp(dir=basedir))
            tx_file = os.path.join(tmpdir, os.path.basename(fname))
            tx_files.append(tx_file)
            orig_files.append(fname)
    return tx_files, orig_files
