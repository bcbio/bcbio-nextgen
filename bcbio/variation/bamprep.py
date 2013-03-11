"""Provide piped, no disk-IO, BAM preparation for variant calling.
Handles independent analysis of chromosome regions, allowing parallel
runs of this step.
"""
import os
import subprocess

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import shared as pipeshared

def _piped_input_cl(data, region, tmp_dir, out_base_file, mark_duplicates):
    """Retrieve the commandline for streaming input into preparation step.
    If marking duplicates, this requires writing an intermediate file since
    MarkDuplicates uses multiple passed on an input.
    """
    broad_runner = broad.runner_from_config(data["config"])
    chrom, start, end = region
    cl = broad_runner.cl_gatk(["-T", "PrintReads",
                               "-L", "%s:%s-%s" % (chrom, start + 1, end),
                               "-R", data["sam_ref"],
                               "-I", data["work_bam"]], tmp_dir)
    if mark_duplicates:
        sel_file = "%s-select%s" % os.path.splitext(out_base_file)
        if not utils.file_exists(sel_file):
            with file_transaction(sel_file) as tx_out_file:
                cl += ["-o", tx_out_file]
                subprocess.check_call(cl)
        dup_metrics = "%s-dup.dup_metrics" % os.path.splitext(out_base_file)[0]
        cl = broad_runner.cl_picard("MarkDuplicates",
                                    [("INPUT", sel_file),
                                     ("OUTPUT", "/dev/stdout"),
                                     ("METRICS_FILE", dup_metrics),
                                     ("PROGRAM_RECORD_ID", "null"),
                                     ("COMPRESSION_LEVEL", "0"),
                                     ("TMP_DIR", tmp_dir)])
    return " ".join(cl)

def _piped_bamprep_region(data, region, out_file, tmp_dir):
    """Do work of preparing BAM input file on the selected region.
    """
    broad_runner = broad.runner_from_config(data["config"])
    algorithm = data["config"]["algorithm"]
    cl = _piped_input_cl(data, region, tmp_dir, out_file, algorithm.get("mark_duplicates", True))
    print cl
    if algorithm.get("recalibrate", True):
        pass
    if algorithm.get("realign", True):
        pass
    with file_transaction(out_file) as tx_out_file:
        subprocess.check_call("{cl} > {tx_out_file}".format(**locals()), shell=True)

def piped_bamprep(data, region=None, out_file=None):
    """Perform full BAM preparation using pipes to avoid intermediate disk IO.

    Handles de-duplication, recalibration and realignment of original BAMs.
    """
    if region[0] == "nochrom":
        prep_bam = pipeshared.write_nochr_reads(data["work_bam"], out_file)
    elif region[0] == "noanalysis":
        prep_bam = pipeshared.write_noanalysis_reads(data["work_bam"], region[1], out_file)
    else:
        if not utils.file_exists(out_file):
            with utils.curdir_tmpdir() as tmp_dir:
                _piped_bamprep_region(data, region, out_file, tmp_dir)
        prep_bam = out_file
    data["work_bam"] = prep_bam
    return [data]
