"""Pipeline utilities to retrieve FASTQ formatted files for processing.
"""
import os
import glob
import subprocess
import contextlib
import collections

import pysam

from bcbio import broad
from bcbio.utils import file_exists, safe_makedir
from bcbio.distributed.transaction import file_transaction

def get_fastq_files(directory, work_dir, item, fc_name, bc_name=None,
                    config=None):
    """Retrieve fastq files for the given lane, ready to process.
    """
    if item.has_key("files") and bc_name is None:
        names = item["files"]
        if isinstance(names, basestring):
            names = [names]
        files = [x if os.path.isabs(x) else os.path.join(directory, x) for x in names]
    else:
        assert fc_name is not None
        lane = item["lane"]
        if bc_name:
            glob_str = "%s_*%s_%s_*_fastq.txt" % (lane, fc_name, bc_name)
        else:
            glob_str = "%s_*%s*_fastq.txt" % (lane, fc_name)
        files = glob.glob(os.path.join(directory, glob_str))
        files.sort()
        if len(files) > 2 or len(files) == 0:
            raise ValueError("Did not find correct files for %s %s %s %s" %
                    (directory, lane, fc_name, files))
    ready_files = []
    for fname in files:
        if fname.endswith(".gz"):
            cl = ["gunzip", fname]
            subprocess.check_call(cl)
            ready_files.append(os.path.splitext(fname)[0])
        elif fname.endswith(".bam"):
            if _pipeline_needs_fastq(config, item):
                ready_files = convert_bam_to_fastq(fname, work_dir, config)
            else:
                ready_files = [fname]
        else:
            assert os.path.exists(fname), fname
            ready_files.append(fname)
    ready_files = [x for x in ready_files if x is not None]
    return ready_files[0], (ready_files[1] if len(ready_files) > 1 else None)

def _pipeline_needs_fastq(config, item):
    do_align = (config["algorithm"].get("aligner", None) and
                item["algorithm"].get("aligner", True))
    has_multiplex = item.get("multiplex") is not None
    do_split = item["algorithm"].get("align_split_size") is not None
    return has_multiplex or (do_align and not do_split)

def convert_bam_to_fastq(in_file, work_dir, config):
    """Convert BAM input file into FASTQ files.
    """
    out_dir = safe_makedir(os.path.join(work_dir, "fastq_convert"))
    out_files = [os.path.join(out_dir, "{0}_{1}.fastq".format(
                 os.path.splitext(os.path.basename(in_file))[0], x))
                 for x in ["1", "2"]]
    if _is_paired(in_file):
        out1, out2 = out_files
    else:
        out1 = out_files[0]
        out2 = None
    if not file_exists(out1):
        broad_runner = broad.runner_from_config(config)
        broad_runner.run_fn("picard_bam_to_fastq", in_file, out1, out2)
    if os.path.getsize(out2) == 0:
        out2 = None
    return [out1, out2]

def _is_paired(bam_file):
    # XXX need development version of pysam for this to work on
    # fastq files without headers (ie. FastqToSam)
    # Instead return true by default and then check after output
    return True
    with contextlib.closing(pysam.Samfile(bam_file, "rb")) as work_bam:
        for read in bam_file:
            return read.is_paired
