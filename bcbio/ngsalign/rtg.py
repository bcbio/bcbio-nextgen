"""Provide indexing and retrieval of files using Real Time Genomics SDF format.

Prepares a sdf representation of reads suitable for indexed retrieval,
normalizing many different input types.

https://github.com/RealTimeGenomics/rtg-tools
"""
import os
import subprocess

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do

def to_sdf(files, data):
    """Convert a fastq or BAM input into a SDF indexed file.
    """
    # BAM
    if len(files) == 1 and files[0].endswith(".bam"):
        qual = []
        format = ["-f", "sam-pe" if bam.is_paired(files[0]) else "sam-se"]
        inputs = [files[0]]
    # fastq
    else:
        qual = ["-q", "illumina" if dd.get_quality_format(data).lower() == "illumina" else "sanger"]
        format = ["-f", "fastq"]
        if len(files) == 2:
            inputs = ["-l", files[0], "-r", files[1]]
        else:
            assert len(files) == 1
            inputs = [files[0]]
    work_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "align_prep"))
    out_file = os.path.join(work_dir,
                            "%s.sdf" % utils.splitext_plus(os.path.basename(os.path.commonprefix(files)))[0])
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = _rtg_cmd(["rtg", "format", "-o", tx_out_file] + format + qual + inputs)
            do.run(cmd, "Format inputs to indexed SDF")
    return out_file

def _rtg_cmd(cmd):
    return "export RTG_JAVA_OPTS='-Xms500m' && export RTG_MEM=2g && " + " ".join(cmd)

def calculate_splits(sdf_file, split_size):
    """Retrieve
    """
    counts = _sdfstats(sdf_file)["counts"]
    splits = []
    cur = 0
    for i in range(counts // split_size + (0 if counts % split_size == 0 else 1)):
        splits.append("%s-%s" % (cur, min(counts, cur + split_size)))
        cur += split_size
    return splits

def _sdfstats(sdf_file):
    cmd = ["rtg", "sdfstats", sdf_file]
    pairs = []
    counts = []
    lengths = []
    for line in subprocess.check_output(_rtg_cmd(cmd), shell=True).split("\n"):
        if line.startswith("Paired arm"):
            pairs.append(line.strip().split()[-1])
        elif line.startswith("Number of sequences"):
            counts.append(int(line.strip().split()[-1]))
        elif line.startswith("Minimum length"):
            lengths.append(int(line.strip().split()[-1]))
    assert len(set(counts)) == 1, counts
    return {"counts": counts[0], "pairs": pairs, "min_size": min(lengths)}

def min_read_size(sdf_file):
    """Retrieve minimum read size from SDF statistics.
    """
    return _sdfstats(sdf_file)["min_size"]

def is_paired(sdf_file):
    """Check if we have paired inputs in a SDF file.
    """
    pairs = _sdfstats(sdf_file)["pairs"]
    return len(set(pairs)) > 1

def to_fastq_apipe_cl(sdf_file, start=None, end=None):
    """Return a command lines to provide streaming fastq input.

    For paired end, returns a forward and reverse command line. For
    single end returns a single command line and None for the pair.
    """
    cmd = ["rtg", "sdf2fastq", "--no-gzip", "-o", "-"]
    if start is not None:
        cmd += ["--start-id=%s" % start]
    if end is not None:
        cmd += ["--end-id=%s" % end]
    if is_paired(sdf_file):
        out = []
        for ext in ["left", "right"]:
            out.append("<(%s)" % _rtg_cmd(cmd + ["-i", os.path.join(sdf_file, ext)]))
        return out
    else:
        cmd += ["-i", sdf_file]
        return ["<(%s)" % _rtg_cmd(cmd), None]