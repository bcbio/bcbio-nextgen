import os
import shutil
import tempfile
from bcbio.utils import chdir, file_exists, safe_makedir, get_in
from bcbio.bam import is_paired
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.distributed.transaction import tx_tmpdir, file_transaction
from bcbio.pipeline import config_utils
from bcbio.log import logger

def run(data):
    """Quantitaive isoforms expression by express"""
    name = dd.get_sample_name(data)
    in_bam = dd.get_transcriptome_bam(data)
    tophat_index = get_in(data, ('genome_resources', 'rnaseq', 'transcriptome_index', 'tophat'))
    if not tophat_index:
        logger.info("Tophat index not found, skipping running eXpress.")
        return None
    tophat_fa = tophat_index.replace("ver", "fa")
    out_dir = os.path.join(dd.get_work_dir(data), "express", name)
    out_file = os.path.join(out_dir, name + ".xprs")
    safe_makedir(out_dir)
    express = config_utils.get_program("express", data['config'])
    if not in_bam:
        logger.info("Transcriptome-mapped BAM file not found, skipping eXpress.")
        return None
    strand = _set_stranded_flag(in_bam, data['config'])
    if not file_exists(out_file):
        with tx_tmpdir() as tmp_dir:
            chdir(tmp_dir)
            ref_transcript = _do_fasta(tophat_fa)
            cmd = ("{express} {strand} {ref_transcript} {in_bam}")
            do.run(cmd.format(**locals()), "Run express", {})
            shutil.move("results.xprs", out_file)
    eff_count_file = _get_column(out_file, out_file.replace(".xprs", "_eff.counts"), 7)
    tpm_file = _get_column(out_file, out_file.replace("xprs", "tpm"), 14)
    fpkm_file = _get_column(out_file, out_file.replace("xprs","fpkm"), 10)
    return (eff_count_file, tpm_file, fpkm_file)

def _do_fasta(fa_file):
    """parse tophat transcriptome.fa to get the correct format for express"""
    out_file = tempfile.NamedTemporaryFile(delete=False, suffix=".fa").name
    transcript = ""
    with open(fa_file) as in_handle:
        with file_transaction(out_file) as tx_out_file:
            with open(tx_out_file, 'w') as out_handle:
                for line in in_handle:
                    if line.startswith(">"):
                        if transcript == line.split(" ")[1]:
                            raise ValueError("You need to update the Tophat index with "
                                             "bcbio_setup_genome.py to run eXpress. This "
                                             "version has some gene annotations that give "
                                             "eXpress issues.")
                        transcript = line.split(" ")[1]
                        out_handle.write(">%s\n" % transcript)
                    else:
                        out_handle.write(line)
    return out_file

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

def _set_stranded_flag(bam_file, config):
    strand_flag = {"unstranded": "",
                   "firststrand": "--rf-stranded",
                   "secondstrand": "--fr-stranded",
                   "firststrand-s": "--r-stranded",
                   "secondstrand-s": "--f-stranded"}
    stranded = get_in(config, ("algorithm", "strandedness"), "unstranded").lower()
    assert stranded in strand_flag, ("%s is not a valid strandedness value. "
            "Valid values are 'firststrand', "
            "'secondstrand' and 'unstranded" % (stranded))
    if not is_paired(bam_file):
        stranded += "-s"
    flag = strand_flag[stranded]
    return flag
