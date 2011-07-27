"""Next gen sequence alignments with Bowtie (http://bowtie-bio.sourceforge.net).
"""
import os
import subprocess

from bcbio.utils import file_transaction

galaxy_location_file = "bowtie_indices.loc"

def align(fastq_file, pair_file, ref_file, out_base, align_dir, config):
    """Before a standard or paired end alignment with bowtie.
    """
    out_file = os.path.join(align_dir, "%s.sam" % out_base)
    if not os.path.exists(out_file):
        cl = [config["program"]["bowtie"]]
        cl += ["-q", "--solexa1.3-quals",
               "-v", config["algorithm"]["max_errors"],
               "-k", 1,
               "-X", 1000, # matches bwa sampe default size
               "-M", 1,
               "--best",
               "--strata",
               "--sam",
               ref_file]
        if pair_file:
            cl += ["-1", fastq_file, "-2", pair_file]
        else:
            cl += [fastq_file]
        cl += [out_file]
        cl = [str(i) for i in cl]
        with file_transaction(out_file):
            subprocess.check_call(cl)
    return out_file

