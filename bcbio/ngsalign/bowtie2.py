"""Next gen sequence alignments with Bowtie2.

http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
"""
import os

from bcbio.pipeline import config_utils
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do

def _bowtie2_args_from_config(config):
    """Configurable high level options for bowtie2.
    """
    qual_format = config["algorithm"].get("quality_format", "")
    if qual_format.lower() == "illumina":
        qual_flags = ["--phred64-quals"]
    else:
        qual_flags = []
    num_cores = config["algorithm"].get("num_cores", 1)
    core_flags = ["-p", str(num_cores)] if num_cores > 1 else []
    return core_flags + qual_flags

def align(fastq_file, pair_file, ref_file, out_base, align_dir, config,
          extra_args=None, names=None):
    """Alignment with bowtie2.
    """
    out_file = os.path.join(align_dir, "%s.sam" % out_base)
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            cl = [config_utils.get_program("bowtie2", config)]
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
            do.run(cl, "Aligning %s and %s with Bowtie2." % (fastq_file, pair_file),
                   None)
    return out_file

# Optional galaxy location file. Falls back on remap_index_fn if not found
galaxy_location_file = "bowtie2_indices.loc"

def remap_index_fn(ref_file):
    """Map sequence references to equivalent bowtie2 indexes.
    """
    return os.path.splitext(ref_file)[0].replace("/seq/", "/bowtie2/")
