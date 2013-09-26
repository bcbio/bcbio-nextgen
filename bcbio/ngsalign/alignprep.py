"""Prepare read inputs (fastq, gzipped fastq and BAM) for parallel NGS alignment.
"""
import copy
import os
import subprocess

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.provenance import do

def create_inputs(data):
    """Index input reads and prepare groups of reads to process concurrently.

    Allows parallelization of alignment beyond processors available on a single
    machine. Uses gbzip and grabix to prepare an indexed fastq file.
    """
    # skip skipping on samples without input files
    if data["files"][0] is None or data["algorithm"].get("align_split_size") is None:
        return [[data]]
    ready_files = _prep_grabix_indexes(data["files"], data["dirs"], data["config"])
    data["files"] = ready_files
    splits = _find_read_splits(ready_files[0], data["algorithm"]["align_split_size"])
    if len(splits) == 1:
        return [[data]]
    else:
        out = []
        for split in splits:
            cur_data = copy.deepcopy(data)
            cur_data["align_split"] = split
            out.append([cur_data])
        return out

# ## merge

def merge_split_alignments(samples, run_parallel):
    """Manage merging split alignments back into a final working BAM file.
    """
    ready = []
    for data in (xs[0] for xs in samples):
        if data.get("align_split"):
            raise NotImplementedError
        else:
            ready.append([data])
    return ready

# ## determine file sections

def _find_read_splits(in_file, split_size):
    """Determine sections of fastq files to process in splits.

    Assumes a 4 line order to input files (name, read, name, quality).
    grabix is 1-based inclusive, so return coordinates in that format.
    """
    gbi_file = in_file + ".gbi"
    with open(gbi_file) as in_handle:
        in_handle.next() # throw away
        num_lines = int(in_handle.next().strip())
    assert num_lines % 4 == 0, "Expected lines to be multiple of 4"
    split_lines = split_size * 4
    chunks = []
    for chunki in range(num_lines // split_lines + min(1, num_lines % split_lines)):
        chunks.append((chunki * split_lines + 1, min((chunki + 1) * split_lines, num_lines)))
    return chunks

# ## bgzip and grabix

def _prep_grabix_indexes(in_files, dirs, config):
    if in_files[0].endswith(".bam") and in_files[1] is None:
        raise NotImplementedError("Prepare BAM indexed files.")
    else:
        return [_grabix_from_fastq(x, dirs, config) if x else None for x in in_files]

def _grabix_from_fastq(in_file, dirs, config):
    """Prepare a bgzipped grabix indexed file from a fastq input.
    """
    grabix = config_utils.get_program("grabix", config)
    if in_file.endswith(".gz"):
        needs_bgzip, needs_gunzip = _check_gzipped_input(in_file, grabix)
    else:
        needs_bgzip, needs_gunzip = True, False
    if needs_bgzip or needs_gunzip:
        out_file = _bgzip_file(in_file, dirs, config, needs_bgzip, needs_gunzip)
    else:
        out_file = in_file
    gbi_file = out_file + ".gbi"
    if not utils.file_exists(gbi_file):
        do.run([grabix, "index", out_file], "Index input with grabix")
    return out_file

def _bgzip_file(in_file, dirs, config, needs_bgzip, needs_gunzip):
    """Handle bgzip of input file, potentially gunzipping an existing file.
    """
    work_dir = utils.safe_makedir(os.path.join(dirs["work"], "align_prep"))
    out_file = os.path.join(work_dir, os.path.basename(in_file) +
                            (".gz" if not in_file.endswith(".gz") else ""))
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            bgzip = config_utils.get_program("bgzip", config)
            assert needs_bgzip
            if needs_gunzip:
                gunzip_cmd = "gunzip -c {in_file} |".format(**locals())
                bgzip_in = "/dev/stdin"
            else:
                gunzip_cmd = ""
                bgzip_in = in_file
            do.run("{gunzip_cmd} {bgzip} -c {bgzip_in} > {tx_out_file}".format(**locals()),
                   "bgzip input file")
    return out_file

def _check_gzipped_input(in_file, grabix):
    """Determine if a gzipped input file is blocked gzip or standard.
    """
    is_bgzip = subprocess.check_output([grabix, "check", in_file])
    if is_bgzip.strip() == "yes":
        return False, False
    else:
        return True, True
