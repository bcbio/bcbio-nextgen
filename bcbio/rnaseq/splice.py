"""High level parallel alternative splice callers.
"""
import os
import copy
import commands
from bcbio.log import logger
from bcbio import bam, utils
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import fastq
from bcbio.rnaseq import rmats
# from bcbio.pipeline import region


def get_callers():
    from bcbio.rnaseq import rmats
    return {"rmats": rmats.run}

def splicecall_prepare(data, run_parallel):
    """Entry point for doing alternative splice callers"""
    gtf_file = dd.get_gtf_file(data)
    caller_fns = get_callers()
    to_process = []
    caller = "rmats"
    for sample in data:
        if dd.get_replicate(sample[0]) == 1:
            mimic = copy.copy(sample[0])
            if caller in dd.get_splicecaller(sample[0]):
                if caller in caller_fns and dd.get_phenotype(mimic) != "control":
                    mimic["rmats_fn"] = caller
                    name = dd.get_sample_name(mimic)
                    rep_mimic = _get_replicate_samples(mimic, data)
                    mimic = _get_paired_samples(mimic, data)
                    if mimic:
                        to_process.append(mimic)
                    else:
                        logger.info("Skipping alternative splice calling. No input sample for %s" % name)
    if to_process:
        after_process = run_parallel("splicecalling", to_process)
        data = _sync(data, after_process)
    return data

def calling(data):
    """Main function to parallelize peak calling."""
    chip_bam = dd.get_work_bam(data)
    if data["work_bam_rep"] != "":
        rep_bam = data["work_bam_rep"]
    else:
        rep_bam = ""
    input_bam = data["work_bam_input"]
    caller_fn = get_callers()[data["rmats_fn"]]
    name = dd.get_sample_name(data)
    myCmd = 'samtools view '+chip_bam+' | head -n 1'
    status,output=commands.getstatusoutput(myCmd)
    read_len=len(output.strip().split('\t')[9])
    fastq_file = fastq.get_fastq_files(data)
    if len(fastq_file) > 1:
        read_pair = "paired"
    else:
        read_pair = "single"
    out_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), data["rmats_fn"], name ))
    out_file = caller_fn(name, chip_bam, rep_bam, input_bam, dd.get_gtf_file(data), out_dir, read_len, read_pair, data["config"])
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
    input_bam = ""
    for origin in data:
        if  dd.get_phenotype(origin[0]) == "control":
            if input_bam != "":
                input_bam = input_bam + "," + dd.get_work_bam(origin[0])
            else:
                input_bam = dd.get_work_bam(origin[0])
    sample["work_bam_input"] = input_bam
    return [sample]

def _get_replicate_samples(sample, data):
    """Get input sample for each chip bam file."""
    dd.get_phenotype(sample)
    rep_bam = ""
    for origin in data:
        if  dd.get_batch(sample) in dd.get_batch(origin[0]) and dd.get_phenotype(sample) in dd.get_phenotype(origin[0]) and dd.get_work_bam(sample) != dd.get_work_bam(origin[0]) and dd.get_phenotype(origin[0]) != "control":
            if rep_bam != "":
                rep_bam = rep_bam + "," + dd.get_work_bam(origin[0])
            else:
                rep_bam = dd.get_work_bam(origin[0])
    sample["work_bam_rep"] = dd.get_work_bam(origin[0])
    return [sample]

def _get_multiplier(samples):
    """Get multiplier to get jobs
       only for samples that have input
    """
    to_process = 1.0
    for sample in samples:
        if dd.get_phenotype(sample[0]) != "control" and dd.get_replicate(sample[0]) == 1:
            to_process += 1.0
    if to_process / len(samples) < 1.0:
        to_process = 1.0
    else:
        to_process = to_process / len(samples)
    return to_process
