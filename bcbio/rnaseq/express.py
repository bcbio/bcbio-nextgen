import os
import shutil
import tempfile
from bcbio.utils import chdir, file_exists, safe_makedir, get_in
from bcbio.bam import is_paired
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.distributed.transaction import tx_tmpdir, file_transaction
from bcbio.pipeline import config_utils
from bcbio.rnaseq import gtf
from bcbio.log import logger

def run(data):
    """Quantitaive isoforms expression by eXpress"""
    name = dd.get_sample_name(data)
    in_bam = dd.get_transcriptome_bam(data)
    if not in_bam:
        logger.info("Transcriptome-mapped BAM file not found, skipping eXpress.")
        return None
    gtf_fasta = gtf.gtf_to_fasta(dd.get_gtf_file(data), dd.get_ref_file(data))
    out_dir = os.path.join(dd.get_work_dir(data), "express", name)
    out_file = os.path.join(out_dir, name + ".xprs")
    safe_makedir(out_dir)
    express = config_utils.get_program("express", data['config'])
    strand = _set_stranded_flag(in_bam, data)
    if not file_exists(out_file):
        with tx_tmpdir() as tmp_dir:
            chdir(tmp_dir)
            cmd = ("{express} {strand} {gtf_fasta} {in_bam}")
            do.run(cmd.format(**locals()), "Run express", {})
            shutil.move("results.xprs", out_file)
    eff_count_file = _get_column(out_file, out_file.replace(".xprs", "_eff.counts"), 7)
    tpm_file = _get_column(out_file, out_file.replace("xprs", "tpm"), 14)
    fpkm_file = _get_column(out_file, out_file.replace("xprs","fpkm"), 10)
    return (eff_count_file, tpm_file, fpkm_file)

def _get_column(in_file, out_file, column):
    """Subset one column from a file
    """
    with file_transaction(out_file) as tx_out_file:
        with open(in_file) as in_handle:
            with open(tx_out_file, 'w') as out_handle:
                for line in in_handle:
                    cols = line.strip().split("\t")
                    if line.find("eff_count") > 0:
                        continue
                    number = cols[column]
                    if column == 7:
                        number = int(round(float(number), 0))
                    out_handle.write("%s\t%s\n" % (cols[1], number))
    return out_file

def _set_stranded_flag(bam_file, data):
    strand_flag = {"unstranded": "",
                   "firststrand": "--rf-stranded",
                   "secondstrand": "--fr-stranded",
                   "firststrand-s": "--r-stranded",
                   "secondstrand-s": "--f-stranded"}
    stranded = dd.get_strandedness(data)
    assert stranded in strand_flag, ("%s is not a valid strandedness value. "
            "Valid values are 'firststrand', "
            "'secondstrand' and 'unstranded" % (stranded))
    if stranded != "unstranded" and not is_paired(bam_file):
        stranded += "-s"
    flag = strand_flag[stranded]
    return flag
