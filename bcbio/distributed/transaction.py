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
import time

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
        config_tmpdir = utils.safe_makedir(os.path.expandvars(config_tmpdir))
        config_tmpdir = os.path.normpath(os.path.join(os.getcwd(), config_tmpdir))
        tmp_dir_base = os.path.join(config_tmpdir, "bcbiotx", str(uuid.uuid4()))
        unique_attempts = 0
        while os.path.exists(tmp_dir_base):
            if unique_attempts > 5:
                break
            tmp_dir_base = os.path.join(config_tmpdir, "bcbiotx", str(uuid.uuid4()))
            time.sleep(1)
            unique_attempts += 1
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
            for dname in [tmp_dir, tmp_dir_base if config_tmpdir else None]:
                if dname and os.path.exists(dname):
                    try:
                        shutil.rmtree(dname, ignore_errors=True)
                    except:
                        pass

@contextlib.contextmanager
def file_transaction(*data_and_files):
    """Wrap file generation in a transaction, moving to output if finishes.

    The initial argument can be the world descriptive `data` dictionary, or
    a `config` dictionary. This is used to identify global settings for
    temporary directories to create transactional files in.
    """
    exts = {".vcf": ".idx", ".bam": ".bai", ".vcf.gz": ".tbi", ".bed.gz": ".tbi"}
    with _flatten_plus_safe(data_and_files) as (safe_names, orig_names):
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
                    utils.safe_makedir(os.path.dirname(orig))
                    # if we are rolling back a directory and it already exists
                    # this will avoid making a nested set of directories
                    if os.path.isdir(orig) and os.path.isdir(safe):
                        shutil.rmtree(orig)
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

@contextlib.contextmanager
def _flatten_plus_safe(data_and_files):
    """Flatten names of files and create temporary file names.
    """
    data_and_files = [x for x in data_and_files if x]
    if isinstance(data_and_files[0], dict):
        data = data_and_files[0]
        rollback_files = data_and_files[1:]
    else:
        data = None
        rollback_files = data_and_files
    tx_files, orig_files = [], []
    base_fname = rollback_files[0]
    if isinstance(base_fname, (list, tuple)):
        base_fname = base_fname[0]
    with tx_tmpdir(data, os.path.dirname(base_fname)) as tmpdir:
        for fnames in rollback_files:
            if isinstance(fnames, basestring):
                fnames = [fnames]
            for fname in fnames:
                tx_file = os.path.join(tmpdir, os.path.basename(fname))
                tx_files.append(tx_file)
                orig_files.append(fname)
        yield tx_files, orig_files
