"""Pipeline utilities to retrieve FASTQ formatted files for processing.
"""
import os
import glob
import subprocess
import contextlib

import pysam

from bcbio import broad
from bcbio.bam import cram
from bcbio.pipeline import alignment
from bcbio.utils import file_exists, safe_makedir

def needs_fastq_conversion(item, config):
    """Check if an item needs conversion to fastq files.
    """
    files = item.get("files", [])
    if isinstance(files, basestring):
        files = [files]
    for f in files:
        if f.endswith(".bam") and _pipeline_needs_fastq(config, item):
            return True
    return False

def get_fastq_files(item):
    """Retrieve fastq files for the given lane, ready to process.
    """
    fastq_dir = item["dirs"]["fastq"]
    if "files" in item:
        names = item["files"]
        if isinstance(names, basestring):
            names = [names]
        files = [x if os.path.isabs(x) else os.path.join(fastq_dir, x) for x in names]
    elif "vrn_file" in item:
        files = []
    else:
        assert item["upload"].get("fc_name") is not None
        lane = item["lane"]
        glob_str = "%s_*%s*_fastq.txt" % (lane, item["upload"]["fc_name"])
        files = glob.glob(os.path.join(fastq_dir, glob_str))
        files.sort()
        if len(files) > 2 or len(files) == 0:
            raise ValueError("Did not find correct files for %s %s %s %s" %
                             (fastq_dir, lane, item["upload"]["fc_name"], files))
    ready_files = []
    for fname in files:
        if fname.endswith(".gz") and _pipeline_needs_fastq(item["config"], item):
            cl = ["gunzip", fname]
            subprocess.check_call(cl)
            ready_files.append(os.path.splitext(fname)[0])
        elif fname.endswith(".bam"):
            if _pipeline_needs_fastq(item["config"], item):
                ready_files = convert_bam_to_fastq(fname, item["dirs"]["work"],
                                                   item, item["dirs"], item["config"])
            else:
                ready_files = [fname]
        else:
            assert os.path.exists(fname), fname
            ready_files.append(fname)
    ready_files = [x for x in ready_files if x is not None]
    return ((ready_files[0] if len(ready_files) > 0 else None),
            (ready_files[1] if len(ready_files) > 1 else None))

def _pipeline_needs_fastq(config, item):
    """Determine if the pipeline can proceed with a BAM file, or needs fastq conversion.
    """
    aligner = config["algorithm"].get("aligner")
    has_multiplex = item.get("multiplex") is not None
    do_split = config["algorithm"].get("align_split_size") is not None
    support_bam = aligner in alignment.metadata.get("support_bam", [])
    return (has_multiplex or
            (aligner and not do_split and not support_bam))

def convert_bam_to_fastq(in_file, work_dir, item, dirs, config):
    """Convert BAM input file into FASTQ files.
    """
    out_dir = safe_makedir(os.path.join(work_dir, "fastq_convert"))

    qual_bin_method = config["algorithm"].get("quality_bin")
    if (qual_bin_method == "prealignment" or
         (isinstance(qual_bin_method, list) and "prealignment" in qual_bin_method)):
        _, sam_ref = alignment.get_genome_ref(item["genome_build"], None, dirs["galaxy"])
        out_bindir = safe_makedir(os.path.join(out_dir, "qualbin"))
        in_file = cram.illumina_qual_bin(in_file, sam_ref, out_bindir, config)

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
        for read in work_bam:
            return read.is_paired
