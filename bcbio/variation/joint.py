"""Joint variant calling with multiple samples: aka squaring off, or backfilling.

Handles the N+1 problem of variant calling by combining and recalling samples
previously calling individually (or in smaller batches). Recalls at all positions found
variable in any of the input samples within each batch. Takes a general approach supporting
GATK's incremental joint discovery (http://www.broadinstitute.org/gatk/guide/article?id=3893)
and FreeBayes's N+1 approach (https://groups.google.com/d/msg/freebayes/-GK4zI6NsYY/Wpcp8nt_PVMJ)
as implemented in bcbio.variation.recall (https://github.com/chapmanb/bcbio.variation.recall).
"""
import toolz as tz

from bcbio.distributed.split import grouped_parallel_split_combine
from bcbio.variation import multi

def _split_by_callable_region(data):
    """Split by callable or variant regions.

    We expect joint calling to be deep in numbers of samples per region, so prefer
    splitting aggressively.
    """
    print tz.get_in(("config", "algorithm", "callable_regions"), data)
    print tz.get_in(("config", "algorithm", "variant_regions"), data)

def square_off(samples, run_parallel):
    """Perform joint calling at all variants within a batch.
    """
    # XXX work in progress
    return samples

    to_process = []
    extras = []
    for data in [x[0] for x in samples]:
        jointcaller = tz.get_in(("config", "algorithm", "jointcaller"), data)
        batch = tz.get_in(("metadata", "batch"), data)
        if jointcaller and batch:
            to_process.append(data)
        else:
            extras.append([data])

    processed = grouped_parallel_split_combine(to_process, _split_by_callable_region,
                                               multi.group_batches_joint, run_parallel,
                                               "square_batch_region", "concat_variant_files",
                                               "vrn_file", ["region", "sam_ref", "config"])
    return processed + extras

def square_batch_region(data, region, align_bams, vrn_files, out_file):
    """Perform squaring of a batch in a supplied region, with input BAMs
    """
    pass
