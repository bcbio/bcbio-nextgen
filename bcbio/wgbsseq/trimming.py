import os

from bcbio.utils import append_stem, replace_directory
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
from bcbio import utils
from bcbio.pipeline import datadict as dd
from bcbio.qc import fastqc
from bcbio.bam import fastq
from bcbio.pipeline import config_utils


def trim(data):
    """Remove adapter for bisulphite conversion sequencing data"""
    in_files = data["files"]
    names = dd.get_sample_name(data)
    work_dir = os.path.join(dd.get_work_dir(data), "trimmed", names)
    out_dir = utils.safe_makedir(work_dir)
    out_files = [
        os.path.join(out_dir,
                     utils.splitext_plus(os.path.basename(in_files[0]))[0]
                     + '_val_1.fq.gz'),
        os.path.join(out_dir,
                     utils.splitext_plus(os.path.basename(in_files[1]))[0]
                     + '_val_2.fq.gz')
    ]

    if utils.file_exists(out_files[0]):
        data["files"] = out_files
        return [[data]]

    trim_galore = config_utils.get_program("trim_galore", data["config"])
    cmd = "{trim_galore} --clip_r1 8 --clip_r2 8 --three_prime_clip_r1 8 --three_prime_clip_r2 8 --length 30 --quality 30 --fastqc --paired -o {tx_out_dir} {files}"
    log_file = os.path.join(out_dir, names + "_cutadapt_log.txt")

    if not utils.file_exists(out_files[0]):
        with file_transaction(out_dir) as tx_out_dir:
            files = "%s %s" % (in_files[0], in_files[1])
            do.run(cmd.format(**locals()), "remove adapters with trimgalore")

    data["files"] = out_files
    return [[data]]


def _run_qc_fastqc(in_files, data, out_dir):
    in_files = fastq.downsample(in_files[0], in_files[1], N=5000000)
    for fastq_file in in_files:
        if fastq_file:
            fastqc.run(fastq_file, data, os.path.join(out_dir, utils.splitext_plus(os.path.basename(fastq_file))[0]))


def _fix_output(in_file, stem, out_dir):
    out_file = utils.splitext_plus(replace_directory(append_stem(in_file, stem), out_dir))
    return  "%s%s" % (out_file[0], out_file[1].replace("fastq", "fq"))
