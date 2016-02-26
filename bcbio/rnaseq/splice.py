"""High level parallel alternative splice callers.
"""
import os
import copy

from bcbio.log import logger
from bcbio import bam, utils
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import fastq
from bcbio.rnaseq import rmats
# from bcbio.pipeline import region


def get_callers():
    from bcbio.rnaseq import rmats
    return {"python RNASeq-MATS.py": rmats.run}

def peakcall_prepare(data, run_parallel):
    """Entry point for doing alternative splice callers"""
    gtf_file = dd.get_gtf_file(data)
    caller_fns = get_callers()
    to_process = []
    for sample in data:
    	fastq_file = get_fastq_files(sample)
    	read_len = bam.fastq.estimate_read_length(fastq_file)
        mimic = copy.copy(sample[0])
        for caller in dd.get_peakcaller(sample[0]):
            if caller in caller_fns and dd.get_phenotype(mimic) == "mutant":
                mimic["rmats_fn"] = caller
                name = dd.get_sample_name(mimic)
                mimic = _get_paired_samples(mimic, data)
                if mimic:
                    to_process.append(mimic)
                else:
                    logger.info("Skipping peak calling. No input sample for %s" % name)
    if to_process:
        after_process = run_parallel("peakcalling", to_process)
        data = _sync(data, after_process)
    return data

def calling(data):
    """Main function to parallelize peak calling."""
    chip_bam = dd.get_work_bam(data)
    input_bam = data["work_bam_input"]
    caller_fn = get_callers()[data["rmats_fn"]]
    name = dd.get_sample_name(data)
    out_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), data["rmats_fn"], name ))
    out_file = caller_fn(name, chip_bam, input_bam, dd.get_genome_build(data), out_dir, data["config"])
    data["rmats_file"] = out_file
    return [[data]]

def _sync(original, processed):
    """
    Add output to data if run sucessfully.
    For now only rmats is available, so no need
    to consider multiple callers.
    """
    for original_sample in original:
        original_sample[0]["rmats_file"] = []
        for processs_sample in processed:
            if dd.get_sample_name(original_sample[0]) == dd.get_sample_name(processs_sample[0]):
                if utils.file_exists(processs_sample[0]["rmats_file"]):
                    original_sample[0]["rmats_file"].append(processs_sample[0]["rmats_file"])
    return original

def _get_paired_samples(sample, data):
    """Get input sample for each chip bam file."""
    dd.get_phenotype(sample)
    for origin in data:
        if  dd.get_batch(sample) in dd.get_batch(origin[0]) and dd.get_phenotype(origin[0]) == "control":
            sample["work_bam_input"] = dd.get_work_bam(origin[0])
            return [sample]

def _get_multiplier(samples):
    """Get multiplier to get jobs
       only for samples that have input
    """
    to_process = 1.0
    for sample in samples:
        if dd.get_phenotype(sample[0]) == "mutant":
            to_process += 1.0
    return to_process / len(samples)
