"""Provide analysis of input files by chromosomal regions.

Handle splitting and analysis of files from chromosomal subsets separated by
no-read regions.
"""
import os

from bcbio.distributed.split import parallel_split_combine
from bcbio import utils

def _split_by_regions(dirname, out_ext, in_key):
    """Split a BAM file data analysis into chromosomal regions.
    """
    def _do_work(data):
        bam_file = data[in_key]
        part_info = []
        nowork = [["nochrom"], ["noanalysis", data["regions"]["noanalysis"]]]
        for region in data["regions"]["analysis"] + nowork:
            out_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], dirname,
                                                      data["name"][-1], region[0]))
            region_outfile = os.path.join(out_dir, "%s%s" %
                                          (os.path.splitext(os.path.basename(bam_file))[0],
                                           out_ext))
            part_info.append((region, region_outfile))
        return None, part_info
    return _do_work

def parallel_prep_region(samples, run_parallel):
    """Perform full pre-variant calling BAM prep work on regions.
    """
    file_key = "work_bam"
    split_fn = _split_by_regions("bamprep", "-prep.bam", file_key)
    return parallel_split_combine(samples, split_fn, run_parallel,
                                  "piped_bamprep", None, file_key, ["config"])

def parallel_variantcall_region(samples, run_parallel):
    """Perform variant calling and post-analysis on samples by region.
    """
    raise NotImplementedError
