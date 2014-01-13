"""Next gen sequence alignments with Bowtie2.

http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
"""
import os
from itertools import ifilter, imap
import pysam
import sys

from bcbio.pipeline import config_utils
from bcbio.utils import file_exists, compose
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio import bam



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

def align(fastq_file, pair_file, ref_file, names, align_dir, data,
          extra_args=None):
    """Alignment with bowtie2.
    """
    config = data["config"]
    analysis_config = ANALYSIS.get(data["analysis"])
    assert analysis_config, "Analysis %s is not supported by bowtie2" % (data["analysis"])
    out_file = os.path.join(align_dir, "%s.sam" % names["lane"])
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            cl = [config_utils.get_program("bowtie2", config)]
            cl += _bowtie2_args_from_config(config)
            cl += extra_args if extra_args is not None else []
            cl += ["-q",
                   "-x", ref_file]
            cl += analysis_config.get("params", [])
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


def filter_multimappers(align_file):
    """
    It does not seem like bowtie2 has a corollary to the -m 1 flag in bowtie,
    there are some options that are close but don't do the same thing. Bowtie2
    sets the XS flag for reads mapping in more than one place, so we can just
    filter on that. This will not work for other aligners.
    """
    type_flag = "b" if bam.is_bam(align_file) else ""
    base, ext = os.path.splitext(align_file)
    align_handle = pysam.Samfile(align_file, "r" + type_flag)
    tmp_out_file = os.path.splitext(align_file)[0] + ".tmp"
    def keep_fn(read):
        return _is_properly_mapped(read) and _is_unique(read)
    keep = ifilter(keep_fn, align_handle)
    with pysam.Samfile(tmp_out_file, "w" + type_flag, template=align_handle) as out_handle:
        for read in keep:
            out_handle.write(read)
    align_handle.close()
    out_handle.close()
    os.rename(tmp_out_file, align_file)
    return align_file

def _is_properly_mapped(read):
    if read.is_paired and not read.is_proper_pair:
        return False
    if read.is_unmapped:
        return False
    return True

def _is_unique(read):
    tags = [x[0] for x in read.tags]
    return "XS" not in tags


ANALYSIS = {"chip-seq": {"params": ["-X", 2000]},
            "RNA-seq": {"params": ["--sensitive", "-X", 2000]}}
