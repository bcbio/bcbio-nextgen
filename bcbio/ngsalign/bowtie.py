"""Next gen sequence alignments with Bowtie (http://bowtie-bio.sourceforge.net).
"""
import os
import subprocess

from bcbio.utils import file_transaction

galaxy_location_file = "bowtie_indices.loc"

def align(fastq_file, pair_file, ref_file, out_base, align_dir, config,
          extra_args=None):
    """Before a standard or paired end alignment with bowtie.
    """
    qual_format = config["algorithm"].get("quality_format", None)
    if qual_format is None or qual_format.lower() == "illumina":
        qual_flags = ["--phred64-quals"]
    else:
        qual_flags = []
    multi_mappers = config["algorithm"].get("multiple_mappers", True)
    multi_flags = ["-M", 1] if multi_mappers else ["-m", 1]
    out_file = os.path.join(align_dir, "%s.sam" % out_base)
    if not os.path.exists(out_file):
        cl = [config["program"]["bowtie"]]
        cl += qual_flags
        cl += multi_flags
        cl += extra_args if extra_args is not None else []
        cl += ["-q",
               "-v", config["algorithm"]["max_errors"],
               "-k", 1,
               "-X", 2000, # default is too selective for most data
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

