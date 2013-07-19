"""
count number of reads mapping to features of transcripts

"""
import os

from bcbio.utils import (which, file_exists, get_in, safe_makedir)
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline.shared import configured_ref_file
from bcbio.pipeline.alignment import sam_to_querysort_sam
from bcbio.provenance import do


def htseq_count(data):
    in_file, gtf_file, out_file = _get_files(data)

    if file_exists(out_file):
        return out_file

    htseq = _choose_htseq_count_executable(data)

    with file_transaction(out_file) as tmp_out_file:
        htseq_cmd = ("{htseq} --mode=union --stranded=no --type=exon "
                     "--idattr=gene_id {in_file} {gtf_file} > {tmp_out_file}")

        cmd = htseq_cmd.format(**locals())
        do.run(cmd, "Running htseq-count on %s." % (in_file),
               data["config"], None)

    return out_file


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
    htseq = get_in(data["config"], ("resources", "htseq-count", "cmd"),
                   "htseq-count")
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
