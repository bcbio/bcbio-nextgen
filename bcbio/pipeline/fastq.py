"""Pipeline utilities to retrieve FASTQ formatted files for processing.
"""
import os
import glob
import subprocess

from bcbio import bam, broad
from bcbio.bam import cram
from bcbio.pipeline import alignment
from bcbio.utils import file_exists, safe_makedir, flatten
from bcbio.distributed.transaction import file_transaction

def needs_fastq_conversion(item, config):
    """Check if an item needs conversion to fastq files.
    """
    if item.get("test_run", False):
        return True
    for f in item.get("files", []):
        if f.endswith(".bam") and _pipeline_needs_fastq(config, item):
            return True
    return False

def get_fastq_files(item):
    """Retrieve fastq files for the given lane, ready to process.
    """
    if "files" in item:
        files = item["files"]
    elif "vrn_file" in item:
        files = []
    else:
        assert item["upload"].get("fc_name") is not None
        fastq_dir = item["dirs"]["fastq"]
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
            fastq_dir = os.path.join(item["dirs"]["work"], "fastq")
            safe_makedir(fastq_dir)
            out_file = os.path.join(fastq_dir,
                                    os.path.basename(os.path.splitext(fname)[0]))
            with file_transaction(out_file) as tx_out_file:
                cmd = "gunzip -c {fname} > {tx_out_file}".format(**locals())
                with open(tx_out_file, "w") as out_handle:
                    subprocess.check_call(cmd, shell=True)
            ready_files.append(out_file)
        elif fname.endswith(".bam"):
            if _pipeline_needs_fastq(item["config"], item):
                ready_files = _convert_bam_to_fastq(fname, item["dirs"]["work"],
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
    support_bam = aligner in alignment.metadata.get("support_bam", [])
    return (has_multiplex or
            (aligner and not support_bam))

def _convert_bam_to_fastq(in_file, work_dir, item, dirs, config):
    """Convert BAM input file into FASTQ files.
    """
    out_dir = safe_makedir(os.path.join(work_dir, "fastq_convert"))

    qual_bin_method = config["algorithm"].get("quality_bin")
    if (qual_bin_method == "prealignment" or
         (isinstance(qual_bin_method, list) and "prealignment" in qual_bin_method)):
        out_bindir = safe_makedir(os.path.join(out_dir, "qualbin"))
        in_file = cram.illumina_qual_bin(in_file, item["sam_ref"], out_bindir, config)

    out_files = [os.path.join(out_dir, "{0}_{1}.fastq".format(
                 os.path.splitext(os.path.basename(in_file))[0], x))
                 for x in ["1", "2"]]
    if bam.is_paired(in_file):
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
