import os
import pysam

from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
from bcbio.utils import file_exists, splitext_plus
from bcbio import bam
from bcbio.log import logger

def fix_umi_dragen_bam(data, bam=None):
    """
    fixes the UMI BAM from DRAGEN. Accepts a pre UMI collapsed BAM file and
    adds several missing tags needed to use fgbio's UMI tools
    """
    if not bam:
        bam = dd.get_work_bam(data)
    base, ext = os.path.splitext(bam)
    sample_name = dd.get_sample_name(data)
    out_bam = os.path.join(dd.get_work_dir(data), "align",
                           sample_name, sample_name + "-fgbio" + ext)
    out_bam = add_fgbio_tags(bam, out_bam, data)
    data = dd.set_work_bam(data, out_bam)
    return data

def add_fgbio_tags(in_bam, out_bam, data):
    """
    Add missing MC (mate cigar), MQ (mate quality) and RX (UMI) tags to a DRAGEN BAM file
    """
    if has_fgbio_tags(in_bam):
        logger.info(f"{in_bam} is properly formatted to run fgbio.")
        return in_bam

    if file_exists(out_bam):
        return out_bam

    samtools = config_utils.get_program("samtools", data)
    bamsormadup = config_utils.get_program("bamsormadup", data)

    with file_transaction(out_bam) as tx_out_bam:
        tx_out_dup = "%s-markdup" % splitext_plus(tx_out_bam)[0]
        cmd = (f"{samtools} view -b {in_bam} | "
               f"{bamsormadup} outputformat=sam tmpfile={tx_out_dup} | "
               f"awk '/^@/ {{print;next}} {{N=split($1,n,\":\");print $0 \"\tRX:Z:\" n[N]}}' | "
               f"{samtools} view -bh > {tx_out_bam}")
        message = f"Adding MC/MQ/RX tags to DRAGEN BAM file {in_bam}."
        do.run(cmd, message)
    return out_bam

def has_fgbio_tags(in_bam):
    """
    Checks to see if the MC/MQ or RX tag are present
    """
    FGBIO_TAGS = set(["MC", "MQ", "RX"])
    return bam.has_tags(in_bam, FGBIO_TAGS)
