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

def parallel_split_combine(args, split_fn, parallel_fn,
                           parallel_name, combine_name,
                           file_key, combine_arg_keys):
    """Split, run split items in parallel then combine to output file.

    split_fn: Split an input file into parts for processing. Returns
      the name of the combined output file along with the individual
      split output names and arguments for the parallel function.
    parallel_fn: Reference to run_parallel function that will run
      single core, multicore, or distributed as needed.
    """
    split_args, combine_map, finished_out = _get_split_tasks(args, split_fn, file_key)
    split_output = parallel_fn(parallel_name, split_args)
    combine_args, final_args = _organize_output(split_output, combine_map,
                                                file_key, combine_arg_keys)
    parallel_fn(combine_name, combine_args)
    return finished_out + final_args

def _organize_output(output, combine_map, file_key, combine_arg_keys):
    """Combine output details for parallelization.

    file_key is the key name of the output file used in merging. We extract
    this file from the output data.
    """
    out_map = collections.defaultdict(list)
    extra_args = {}
    final_args = []
    for data in output:
        cur_file = data[file_key]
        cur_out = combine_map[cur_file]
        out_map[cur_out].append(cur_file)
        extra_args[cur_out] = [data[x] for x in combine_arg_keys]
        data[file_key] = cur_out
        if data not in final_args:
            final_args.append([data])
    combine_args = [[v, k] + extra_args[k] for (k, v) in out_map.iteritems()]
    return combine_args, final_args

def _get_split_tasks(args, split_fn, file_key):
    """Split up input files and arguments, returning arguments for parallel processing.
    """
    split_args = []
    combine_map = {}
    finished_out = []
    for data in args:
        out_final, out_parts = split_fn(*data)
        for parts in out_parts:
            split_args.append(data + list(parts))
        for part_file in [x[-1] for x in out_parts]:
            combine_map[part_file] = out_final
        if len(out_parts) == 0:
            data[0][file_key] = out_final
            finished_out.append(data)
    return split_args, combine_map, finished_out
