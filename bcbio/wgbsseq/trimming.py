import os
import sys
import glob
import os.path as op
import shutil
from collections import Counter

from bcbio.utils import (file_exists, append_stem, replace_directory, symlink_plus)
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio import utils
from bcbio.pipeline import datadict as dd
from bcbio.qc import fastqc
from bcbio.bam import fastq
from bcbio.pipeline import config_utils
from bcbio.log import logger


def trimming(data):
    """
    Remove adapter for bisulphite conversion sequencing data
    """
    in_file = data["files"]
    names = dd.get_sample_name(data)
    paired = is_rrbs = is_directional = ""
    if is_rrbs:
        is_rrbs = "--rrbs"
    if not is_directional and is_rrbs:
        is_directional = "--non_directional"

    work_dir = os.path.join(dd.get_work_dir(data), "trimmed", names)
    out_dir = utils.safe_makedir(work_dir)

    _run_qc_fastqc(in_file, data, op.join(out_dir, "before"))

    if len(in_file) == 1:
        out_files = [op.join(out_dir, names + "_1.fastq.gz"), None]
        tmp_files = [op.join(out_dir, names + "_tmp_1.fastq.gz"), None]
    else:
        out_files = [op.join(out_dir, names + "_1.fastq.gz"),
                     op.join(out_dir, names + "_2.fastq.gz")]
        tmp_files = [op.join(out_dir, names + "_tmp_1.fastq.gz"),
                     op.join(out_dir, names + "_tmp_2.fastq.gz")]
    if utils.file_exists(out_files[0]):
        data["files"] = out_files
        return [[data]]
    trim_galore = config_utils.get_program("trim_galore", data["config"])
    cmd = "{trim_galore} {is_directional} {is_rrbs} --length 30 --quality 30 {paired} -o {tx_dir} {files}"
    log_file = op.join(out_dir, names + "_cutadapt_log.txt")
    # cutadapt = config_utils.get_program("cutadapt", data["config"])
    # cmd = "{cutadapt} -q 30 -m 25 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o {first_tmp} -p {second_tmp} {files} | tee > {log_file}"
    if not utils.file_exists(out_files[0]):
        with file_transaction(out_files) as txs:
            files = "%s %s" % (in_file[0], in_file[1])
            do.run(cmd.format(**locals()), "remove adapters with trimgalore.")

    # cmd = "{cutadapt} -m 25 -u -6 -u 8 -U 10 -o {first_out} -p {second_out} {files}"
    # if not utils.file_exists(out_files[0]):
    #     with file_transaction(out_files) as txs:
    #         files = "%s %s" % (tmp_files[0], tmp_files[1])
    #         first_out, second_out = txs
    #         do.run(cmd.format(**locals()), "remove 6 nts, 2 step")
    #     [utils.remove_safe(tmp_file) for tmp_file in tmp_files]
    data["files"] = out_files
    _run_qc_fastqc(out_files, data, op.join(out_dir, "after"))
    return [[data]]

def _run_qc_fastqc(in_files, data, out_dir):
    in_files = fastq.downsample(in_files, None, data=data, N=5000000)
    for fastq_file in in_files:
            if fastq_file:
                fastqc.run(fastq_file, data, op.join(out_dir, utils.splitext_plus(op.basename(fastq_file))[0]), rename=False)

def _fix_output(in_file, stem, out_dir):
    out_file = utils.splitext_plus(replace_directory(append_stem(in_file, stem), out_dir))
    return  "%s%s" % (out_file[0], out_file[1].replace("fastq", "fq"))
