"""Next-gen alignments with BWA (http://bio-bwa.sourceforge.net/)
"""
import os
import subprocess

from bcbio.utils import file_transaction

galaxy_location_file = "bwa_index.loc"

def align(fastq_file, pair_file, ref_file, out_base, align_dir, config):
    """Perform a BWA alignment, generating a SAM file.
    """
    sai1_file = os.path.join(align_dir, "%s_1.sai" % out_base)
    sai2_file = (os.path.join(align_dir, "%s_2.sai" % out_base)
                 if pair_file else None)
    sam_file = os.path.join(align_dir, "%s.sam" % out_base)
    if not os.path.exists(sam_file):
        if not os.path.exists(sai1_file):
            with file_transaction(sai1_file):
                _run_bwa_align(fastq_file, ref_file, sai1_file, config)
        if sai2_file and not os.path.exists(sai2_file):
            with file_transaction(sai2_file):
                _run_bwa_align(pair_file, ref_file, sai2_file, config)
        align_type = "sampe" if sai2_file else "samse"
        sam_cl = [config["program"]["bwa"], align_type, ref_file, sai1_file]
        if sai2_file:
            sam_cl.append(sai2_file)
        sam_cl.append(fastq_file)
        if sai2_file:
            sam_cl.append(pair_file)
        with file_transaction(sam_file):
            with open(sam_file, "w") as out_handle:
                subprocess.check_call(sam_cl, stdout=out_handle)
    return sam_file

def _run_bwa_align(fastq_file, ref_file, out_file, config):
    aln_cl = [config["program"]["bwa"], "aln",
              "-n %s" % config["algorithm"]["max_errors"],
              "-k %s" % config["algorithm"]["max_errors"],
              ref_file, fastq_file]
    with open(out_file, "w") as out_handle:
        subprocess.check_call(aln_cl, stdout=out_handle)

