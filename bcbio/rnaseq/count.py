"""
count number of reads mapping to features of transcripts

"""
from bcbio.log import logger
from bcbio.utils import (which, file_exists, get_in, replace_suffix,
                         safe_makedir)
from bcbio.distributed.transaction import file_transaction
import subprocess
import os
import pysam


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

    data["count_file"] = out_file

    return data

def _get_files(data):
    in_file = _get_sam_file(data)
    config = data["config"]
    gtf_file = get_in(config, ("algorithm", "transcripts"))
    out_file = _get_out_file(in_file, config)
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

def _gtf_exists(config):
    transcript_file = get_in(config, ("algorithm", "transcripts"))
    return file_exists(transcript_file)

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

def is_sam(in_file):
    _, ext = os.path.splitext(in_file)
    if ext == ".sam":
        return True
    else:
        return False

def is_bam(in_file):
    _, ext = os.path.splitext(in_file)
    if ext == ".bam":
        return True
    else:
        return False

def convert_bam_to_sam(in_file):
    if not is_bam(in_file):
        raise ValueError("Non BAM file passed to convert_sam_to_bam: "
                         "%s" % (in_file))
    out_file = replace_suffix(in_file, ".sam")
    if file_exists(out_file):
        return out_file
    with file_transaction(out_file) as tmp_out_file:
        pysam.view("-h", "-o" + tmp_out_file, in_file)
    return out_file

def _get_sam_file(data):
    in_file = data["work_bam"]
    return convert_bam_to_sam(in_file)
