from os import path
import subprocess

from bcbio.pipeline import config_utils
from bcbio.utils import safe_makedir, file_exists
from bcbio.provenance import do

CLEANUP_FILES = ["Aligned.out.sam", "Log.out", "Log.progress.out"]

fastq_file = "/Users/rory/tmp/data/test_fastq_1.fastq"
pair_file = "/Users/rory/tmp/data/test_fastq_2.fastq"
align_dir = "/Users/rory/tmp/star_test"
ref_file = "/Users/rory/tmp/GRCm38/STAR"
out_base = "/Users/rory/tmp/star_test/test"


def align(fastq_file, pair_file, ref_file, out_base, align_dir, config,
          names=None):
    out_prefix = path.join(align_dir, out_base)
    out_file = out_prefix + "Aligned.out.sam"
    if file_exists(out_file):
        return out_file
    star_path = config_utils.get_program("STAR", config)
    fastq = " ".join([fastq_file, pair_file])
    num_cores = config["algorithm"].get("num_cores", 1)
    safe_makedir(align_dir)
    cmd = ("{star_path} --genomeDir {ref_file} --readFilesIn {fastq} "
           "--runThreadN {num_cores} --outFileNamePrefix {out_prefix} "
           "--outReadsUnmapped Fastx --outFilterMultimapNmax 10")
    run_message = "Running STAR aligner on %s and %s." % (pair_file, ref_file)
    do.run(cmd.format(**locals()), run_message, None)
    return out_file


def remap_index_fn(ref_file):
    """Map sequence references to equivalent star indexes
    """
    return path.join(path.dirname(path.dirname(ref_file)), "star")
