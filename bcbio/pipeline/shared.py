"""Pipeline functionality shared amongst multiple analysis types.
"""
import os
import collections
from contextlib import closing

import pysam

from bcbio import broad
from bcbio.pipeline.alignment import get_genome_ref
from bcbio.utils import file_exists, safe_makedir

# ## Split/Combine helpers

def combine_bam(in_files, out_file, config):
    """Parallel target to combine multiple BAM files.
    """
    runner = broad.runner_from_config(config)
    runner.run_fn("picard_merge", in_files, out_file)
    return out_file

def split_bam_by_chromosome(output_ext, file_key):
    """Provide targets to process a BAM file by individual chromosome regions.
    """
    def _do_work(data):
        bam_file = data[file_key]
        out_file = "{base}{ext}".format(base=os.path.splitext(bam_file)[0],
                                        ext=output_ext)
        part_info = []
        if not file_exists(out_file):
            work_dir = safe_makedir(
                "{base}-split".format(base=os.path.splitext(out_file)[0]))
            with closing(pysam.Samfile(bam_file, "rb")) as work_bam:
                for chr_ref in work_bam.references:
                    chr_out = os.path.join(work_dir,
                                           "{base}-{ref}{ext}".format(
                                               base=os.path.splitext(os.path.basename(bam_file))[0],
                                               ref=chr_ref, ext=output_ext))
                    part_info.append((chr_ref, chr_out))
        return out_file, part_info
    return _do_work

# ## Retrieving file information from configuration variables

def configured_ref_file(name, config, sam_ref):
    """Full path to a reference file specified in the configuration.

    Resolves non-absolute paths relative to the base genome reference directory.
    """
    ref_file = config["algorithm"].get(name, None)
    if ref_file:
        if not os.path.isabs(ref_file):
            base_dir = os.path.dirname(os.path.dirname(sam_ref))
            ref_file = os.path.join(base_dir, ref_file)
    return ref_file

def configured_vrn_files(config, sam_ref):
    """Full path to all configured files for variation assessment.
    """
    names = ["dbsnp", "train_hapmap", "train_1000g_omni", "train_indels"]
    VrnFiles = collections.namedtuple("VrnFiles", names)
    return apply(VrnFiles, [configured_ref_file(n, config, sam_ref) for n in names])

def ref_genome_info(info, config, dirs):
    """Retrieve reference genome information from configuration variables.
    """
    genome_build = info.get("genome_build", None)
    (_, sam_ref) = get_genome_ref(genome_build, config["algorithm"]["aligner"],
                                  dirs["galaxy"])
    return genome_build, sam_ref
