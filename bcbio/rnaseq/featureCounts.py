import os

from bcbio.utils import (file_exists, get_in, safe_makedir)
from bcbio.pipeline import config_utils
from bcbio.log import logger
from bcbio.rnaseq.count import htseq_count
from bcbio.bam import is_paired
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction

try:
    import pandas as pd
except ImportError:
    pd = None

def count(data):
    """
    count reads mapping to genes using featureCounts
    falls back on htseq_count method if featureCounts is not
    found
    """
    in_bam = data["work_bam"]
    gtf_file = data["genome_resources"]["rnaseq"]["transcripts"]
    work_dir = data["dirs"].get("work", "work")
    out_dir = os.path.join(work_dir, "htseq-count")
    safe_makedir(out_dir)
    count_file = os.path.join(out_dir, data['rgnames']['sample']) + ".counts"
    if file_exists(count_file):
        return count_file

    config = data["config"]

    try:
        featureCounts = config_utils.get_program("featureCounts", config)
    except config_utils.CmdNotFound:
        logger.info("featureCounts not found, falling back to htseq-count "
                    "for feature counting. You can upgrade the tools to "
                    "install featureCount with bcbio_nextgen.py upgrade "
                    "--tools.")
        return htseq_count(data)

    paired_flag = _paired_flag(in_bam)
    strand_flag = _strand_flag(config)

    cmd = ("{featureCounts} -a {gtf_file} -o {tx_count_file} -s {strand_flag} "
           "{paired_flag} {in_bam}")

    message = ("Count reads in {tx_count_file} mapping to {gtf_file} using "
               "featureCounts")
    with file_transaction(count_file) as tx_count_file:
        do.run(cmd.format(**locals()), message.format(**locals()))
    fixed_count_file = _format_count_file(count_file)
    os.rename(fixed_count_file, count_file)

    return count_file

def _format_count_file(count_file):
    """
    this cuts the count file produced from featureCounts down to
    a two column file of gene ids and number of reads mapping to
    each gene
    """
    COUNT_COLUMN = 5
    out_file = os.path.splitext(count_file)[0] + ".fixed.counts"
    if file_exists(out_file):
        return out_file

    df = pd.io.parsers.read_table(count_file, sep="\t", index_col=0, header=1)
    df_sub = df.ix[:, COUNT_COLUMN]
    with file_transaction(out_file) as tx_out_file:
        df_sub.to_csv(tx_out_file, sep="\t", index_label="id", header=False)
    return out_file


def _strand_flag(config):
    """
    0: unstranded 1: stranded 2: reverse stranded
    """
    strand_flag = {"unstranded": "0",
                   "firststrand": "2",
                   "secondstrand": "1"}
    stranded =  get_in(config, ("algorithm", "strandedness"),
                       "unstranded").lower()

    assert stranded in strand_flag, ("%s is not a valid strandedness value. "
                                     "Valid values are 'firststrand', 'secondstrand', "
                                     "and 'unstranded")
    return strand_flag[stranded]

def _paired_flag(bam_file):
    if is_paired(bam_file):
        return "-p -B -C"
    else:
        return ""



