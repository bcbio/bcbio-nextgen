"""
count number of reads mapping to features of transcripts

"""
from bcbio.log import logger
from bcbio.utils import (which, file_exists, get_in, safe_makedir)
from bcbio.distributed.transaction import file_transaction
import subprocess
import os
from bcbio.pipeline.shared import configured_ref_file
from bcbio.pipeline.alignment import sam_to_querysort_sam
import pandas as pd


def htseq_count(data):
    in_file, gtf_file, out_file = _get_files(data)

    if file_exists(out_file):
        return out_file

    htseq = _choose_htseq_count_executable(data)

    with file_transaction(out_file) as tmp_out_file:
        htseq_cmd = ("{htseq} --mode=union --stranded=no --type=exon "
                     "--idattr=gene_id {in_file} {gtf_file} > {tmp_out_file}")

        cmd = htseq_cmd.format(**locals())
        logger.info("Running htseq-count on {in_file} with command: "
                    "{cmd}".format(**locals()))
        subprocess.check_call(cmd, shell=True)

    return out_file

def combine_htseq_count_files(data):
    count_files = [x.get("count_file", {}) for x in data]
    out_file = os.path.join(os.path.dirname(count_files[0]), "combined.counts")
    if file_exists(out_file):
        return out_file
    df = _combine_htseq_count_files(count_files)
    df.to_csv(out_file, sep="\t", index_label="id")
    return out_file

def _combine_htseq_count_files(files):
    f = files.pop()
    df = pd.io.parsers.read_table(f, sep="\t", index_col=0, header=None,
                                  names=[os.path.basename(f)])
    for f in files:
        df = df.join(pd.io.parsers.read_table(f, sep="\t", index_col=0, header=None,
                                              names=[os.path.basename(f)]))
    return df

def _get_files(data):
    in_file = _get_sam_file(data)
    gtf_file = _get_gtf_file(data)
    out_file = _get_out_file(in_file, data)
    return in_file, gtf_file, out_file

def _get_out_file(in_file, config):
    work_dir = config["dirs"].get("work", "work")
    out_dir = os.path.join(work_dir, "htseq-count")
    safe_makedir(out_dir)
    base, _ = os.path.splitext(os.path.basename(in_file))
    return os.path.join(out_dir, base + ".counts")

def _choose_htseq_count_executable(data):
    htseq = get_in(data["config"], ("resources", "htseq-count", "cmd"), "htseq-count")
    return which(htseq)

def _htseq_is_installed(config):
    if _choose_htseq_count_executable(config):
        return True
    return False

def _get_gtf_file(data):
    ref_file = data["sam_ref"]
    return configured_ref_file("transcripts", data["config"], ref_file)

def _gtf_exists(config):
    gtf_file = _get_gtf_file(config)
    return file_exists(gtf_file)

def is_countfile(in_file):
    with open(in_file) as in_handle:
        firstline = in_handle.next().split("\t")
    if len(firstline) != 2:
        return False
    try:
        int(firstline[1])
    except ValueError:
        return False
    return True


def _get_sam_file(data):
    in_file = data["work_bam"]
    config = data["config"]
    return sam_to_querysort_sam(in_file, config)
