"""Next gen sequence alignments with Bowtie (http://bowtie-bio.sourceforge.net).
"""
import os

from bcbio.pipeline import config_utils
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do

galaxy_location_file = "bowtie_indices.loc"

def _bowtie_args_from_config(config):
    """Configurable high level options for bowtie.
    """
    qual_format = config["algorithm"].get("quality_format", None)
    if qual_format is None or qual_format.lower() == "illumina":
        qual_flags = ["--phred64-quals"]
    else:
        qual_flags = []
    multi_mappers = config["algorithm"].get("multiple_mappers", True)
    multi_flags = ["-M", 1] if multi_mappers else ["-m", 1]
    cores = config.get("resources", {}).get("bowtie", {}).get("cores", None)
    num_cores = config["algorithm"].get("num_cores", 1)
    core_flags = ["-p", str(num_cores)] if num_cores > 1 else []
    return core_flags + qual_flags + multi_flags

def align(fastq_file, pair_file, ref_file, names, align_dir, data,
          extra_args=None):
    """Do standard or paired end alignment with bowtie.
    """
    num_hits = 1
    if data["analysis"].lower().startswith("smallrna-seq"):
        num_hits = 1000
    config = data['config']
    out_file = os.path.join(align_dir, "%s.sam" % names["lane"])
    if not file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cl = [config_utils.get_program("bowtie", config)]
            cl += _bowtie_args_from_config(config)
            cl += extra_args if extra_args is not None else []
            cl += ["-q",
                   "-v", 2,
                   "-k", num_hits,
                   "-X", 2000, # default is too selective for most data
                   "--best",
                   "--strata",
                   "--sam",
                   ref_file]
            if pair_file:
                cl += ["-1", fastq_file, "-2", pair_file]
            else:
                cl += [fastq_file]
            cl += [tx_out_file]
            cl = [str(i) for i in cl]
            do.run(cl, "Running Bowtie on %s and %s." % (fastq_file, pair_file), None)
    return out_file
