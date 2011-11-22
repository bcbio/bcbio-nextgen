"""Next gen sequence alignments with Bowtie2.

http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
"""
import os
import subprocess

from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.ngsalign import bowtie

galaxy_location_file = bowtie.galaxy_location_file

def _bowtie2_args_from_config(config):
    """Configurable high level options for bowtie2.
    """
    qual_format = config["algorithm"].get("quality_format", "")
    if qual_format.lower() == "illumina":
        qual_flags = ["--phred64-quals"]
    else:
        qual_flags = []
    cores = config.get("resources", {}).get("bowtie", {}).get("cores", None)
    core_flags = ["-p", str(cores)] if cores else []
    return core_flags + qual_flags

def align(fastq_file, pair_file, ref_file, out_base, align_dir, config,
          extra_args=None):
    """Alignment with bowtie2.
    """
    out_file = os.path.join(align_dir, "%s.sam" % out_base)
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            cl = [config["program"].get("bowtie2", "bowtie2")]
            cl += _bowtie2_args_from_config(config)
            cl += extra_args if extra_args is not None else []
            cl += ["-q",
                   "--sensitive",
                   "-X", 2000, # default is too selective for most data
                   "-x", ref_file]
            if pair_file:
                cl += ["-1", fastq_file, "-2", pair_file]
            else:
                cl += ["-U", fastq_file]
            cl += ["-S", tx_out_file]
            cl = [str(i) for i in cl]
            subprocess.check_call(cl)
    return out_file

def remap_index_fn(ref_file):
    """Map bowtie references to equivalent bowtie2 indexes.
    """
    return ref_file.replace("/bowtie/", "/bowtie2/")
