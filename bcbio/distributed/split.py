"""Split files or tasks for distributed processing across multiple machines.

This tackles parallel work within the context of a program, where we split
based on input records like fastq or across regions like chromosomes in a
BAM file. Following splitting, individual records and run and then combined
back into a summarized output file.

This provides a framework for that process, making it easier to utilize with
splitting specific code.
"""
import os
import collections

from bcbio.utils import file_transaction, file_exists

def parallel_split_combine(args, split_fn, parallel_fn,
                           parallel_name, combine_name,
                           file_index, extra_combine_args):
    """Split, run split items in parallel then combine to output file.

    split_fn: Split an input file into parts for processing. Returns
      the name of the combined output file along with the individual
      split output names and arguments for the parallel function.
    parallel_fn: Reference to run_parallel function that will run
      single core, multicore, or distributed as needed.
    """
    split_args, combine_map = _get_split_tasks(args, split_fn)
    split_output = parallel_fn(parallel_name, split_args)
    combine_args, final_args = _organize_output(split_output, combine_map,
                                                file_index, extra_combine_args)
    parallel_fn(combine_name, combine_args)
    return final_args

def _organize_output(output, combine_map, file_index, extra_combine_args):
    """Combine output details for parallelization.

    file_index points to the position in the output to find
    the output file. We should be using dictionaries here
    instead.
    """
    out_map = collections.defaultdict(list)
    final_args = []
    for cur in output:
        cur_file = cur[file_index]
        cur_out = combine_map[cur_file]
        out_map[cur_out].append(cur_file)
        cur_arg = list(cur)
        cur_arg[file_index] = cur_out
        if cur_arg not in final_args:
            final_args.append(cur_arg)
    combine_args = [[v, k] + extra_combine_args for (k, v) in out_map.iteritems()]
    return combine_args, final_args

def _get_split_tasks(args, split_fn):
    """Split up input files and arguments, returning arguments for parallel processing.
    """
    split_args = []
    combine_map = {}
    for cur_arg in args:
        out_final, out_parts = split_fn(*cur_arg)
        for parts in out_parts:
            split_args.append(cur_arg + parts)
        for part_file in [x[-1] for x in out_parts]:
            combine_map[part_file] = out_final
    return split_args, combine_map
