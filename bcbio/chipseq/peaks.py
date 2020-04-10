"""High level parallel chip-seq analysis
"""
import os
import copy
import toolz as tz
import subprocess
import shutil
from tempfile import NamedTemporaryFile

import bcbio.bam as bam
from bcbio.log import logger
from bcbio import utils
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
from bcbio.chipseq import atac

def get_callers():
    """Get functions related to each caller"""
    from bcbio.chipseq import macs2
    return {"macs2": macs2.run}

def peakcall_prepare(data, run_parallel):
    """Entry point for doing peak calling"""
    caller_fns = get_callers()
    to_process = []
    for sample in data:
        mimic = copy.copy(sample[0])
        callers = dd.get_peakcaller(sample[0])
        if not isinstance(callers, list):
            callers = [callers]
        for caller in callers:
            if caller in caller_fns:
                mimic["peak_fn"] = caller
                name = dd.get_sample_name(mimic)
                mimic = _check(mimic, data)
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
    method = dd.get_chip_method(data)
    caller_fn = get_callers()[data["peak_fn"]]
    if method == "chip":
        chip_bam = data.get("work_bam")
        input_bam = data.get("work_bam_input", None)
        name = dd.get_sample_name(data)
        out_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), data["peak_fn"], name))
        out_files = caller_fn(name, chip_bam, input_bam, dd.get_genome_build(data), out_dir,
                            dd.get_chip_method(data), data["resources"], data)
        greylistdir = greylisting(data)
        data.update({"peaks_files": out_files})
        if greylistdir:
            data["greylist"] = greylistdir
    if method == "atac":
        fractions = list(atac.ATACRanges.keys()) + ["full"]
        for fraction in fractions:
            MIN_READS_TO_CALL = 1000
            chip_bam = tz.get_in(("atac", "align", fraction), data)
            if bam.count_alignments(chip_bam, data) < MIN_READS_TO_CALL:
                logger.warn(f"{chip_bam} has less than {MIN_READS_TO_CALL}, peak calling will fail so skip this fraction.")
                continue
            logger.info(f"Running peak calling with {data['peak_fn']} on the {fraction} fraction of {chip_bam}.")
            name = dd.get_sample_name(data) + f"-{fraction}"
            out_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), data["peak_fn"], name))
            out_files = caller_fn(name, chip_bam, None, dd.get_genome_build(data), out_dir,
                                  dd.get_chip_method(data), data["resources"], data)
            data = tz.assoc_in(data, ("peaks_files", fraction), out_files)
    return [[data]]

def _sync(original, processed):
    """
    Add output to data if run sucessfully.
    For now only macs2 is available, so no need
    to consider multiple callers.
    """
    for original_sample in original:
        original_sample[0]["peaks_files"] = {}
        for process_sample in processed:
            if dd.get_sample_name(original_sample[0]) == dd.get_sample_name(process_sample[0]):
                for key in ["peaks_files"]:
                    if process_sample[0].get(key):
                        original_sample[0][key] = process_sample[0][key]
    return original

def _check(sample, data):
    """Get input sample for each chip bam file."""
    if dd.get_chip_method(sample).lower() == "atac":
        return [sample]
    if dd.get_phenotype(sample) == "input":
        return None
    for origin in data:
        if dd.get_batch(sample) in (dd.get_batches(origin[0]) or []) and dd.get_phenotype(origin[0]) == "input":
            sample["work_bam_input"] = origin[0].get("work_bam")
            return [sample]
    return [sample]

def _get_multiplier(samples):
    """Get multiplier to get jobs
       only for samples that have input
    """
    to_process = 1.0
    to_skip = 0
    for sample in samples:
        if dd.get_phenotype(sample[0]) == "chip":
            to_process += 1.0
        elif dd.get_chip_method(sample[0]).lower() == "atac":
            to_process += 1.0
        else:
            to_skip += 1.0
    mult = (to_process - to_skip) / len(samples)
    if mult <= 0:
        mult = 1 / len(samples)
    return max(mult, 1)

def greylisting(data):
    """
    Run ChIP-seq greylisting
    """
    input_bam = data.get("work_bam_input", None)
    if not input_bam:
        logger.info("No input control BAM file detected, skipping greylisting.")
        return None
    try:
        greylister = config_utils.get_program("chipseq-greylist", data)
    except config_utils.CmdNotFound:
        logger.info("No greylister found, skipping greylisting.")
        return None
    greylistdir = os.path.join(os.path.dirname(input_bam), "greylist")
    if os.path.exists(greylistdir):
        return greylistdir
    cmd = "{greylister} --outdir {txgreylistdir} {input_bam}"
    message = "Running greylisting on %s." % input_bam
    with file_transaction(greylistdir) as txgreylistdir:
        utils.safe_makedir(txgreylistdir)
        try:
            do.run(cmd.format(**locals()), message)
        except subprocess.CalledProcessError as msg:
            if str(msg).find("Cannot take a larger sample than population when 'replace=False'") >= 0:
                logger.info("Skipping chipseq greylisting because of small sample size: %s"
                            % dd.get_sample_name(data))
                return None
    return greylistdir

def consensus(summitfiles, consensusfile, data, pad=250):
    """call consensus peaks from a set of bed files of summits
    we use this method: https://bedops.readthedocs.io/en/latest/content/usage-examples/master-list.html
    """
    if utils.file_exists(consensusfile):
        return consensusfile

    try:
        bedops = config_utils.get_program("bedops", data)
    except config_utils.CmdNotFound: 
        logger.info("bedops not found, skipping consensus peak calling. do a --tools update to install "
                    "bedops.")
        return None
    try:
        sortbed = config_utils.get_program("sort-bed", data)
    except config_utils.CmdNotFound: 
        logger.info("sort-bed not found, skipping consensus peak calling. do a --tools update to install "
                    "sort-bed.")
        return None
    try:
        bedmap = config_utils.get_program("bedmap", data)
    except config_utils.CmdNotFound: 
        logger.info("bedmap not found, skipping consensus peak calling. do a --tools update to install "
                    "bedmap.")
        return None

    logger.info(f"Calling peaks on {','.join(summitfiles)}")
    with file_transaction(consensusfile) as tx_consensus_file:
        message = f"Combining summits of {' '.join(summitfiles)} and expanding {pad} bases."
        with utils.tmpfile(suffix=".bed") as tmpbed:
            slopcommand = f"{bedops} --range {pad} -u {' '.join(summitfiles)} > {tmpbed}"
            do.run(slopcommand, message)
            iteration = 0
            solutions = []
            while os.path.getsize(tmpbed):
                iteration = iteration + 1
                iterationbed = NamedTemporaryFile(suffix=".bed", delete=False).name
                with utils.tmpfile(suffix="bed") as mergedbed, \
                     utils.tmpfile(suffix="bed") as intermediatebed, \
                     utils.tmpfile(suffix="bed") as leftoverbed:
                    mergecmd = (f"{bedops} -m --range 0:-1 {tmpbed} | "
                                f"{bedops} -u --range 0:1 - > "
                                f"{mergedbed}")
                    message = f"Merging non-overlapping peaks, iteration {iteration}."
                    do.run(mergecmd, message)
                    nitems = len(open(mergedbed).readlines())
                    message = f"Considering {nitems} peaks, choosing the highest score for overlapping peaks."
                    highscorecmd = (f"{bedmap} --max-element {mergedbed} {tmpbed} |"
                                    f"{sortbed} - > "
                                    f"{iterationbed}")
                    do.run(highscorecmd, message)
                    message = f"Checking if there are peaks left to merge."
                    anyleftcmd = (f"{bedops} -n 1 {tmpbed} {iterationbed} > {intermediatebed}")
                    do.run(anyleftcmd, message)
                    shutil.move(intermediatebed, tmpbed)
                    solutions.append(iterationbed)
        message = f"Creating final consensus peak file: {consensusfile}."
        consensuscmd = (f"{bedops} -u {' '.join(solutions)} > {tx_consensus_file}")
        do.run(consensuscmd, message)
    return consensusfile

def call_consensus(samples):
    """call consensus peaks on the summit files from a set of ChiP/ATAC samples
    """
    data = samples[0][0]
    new_samples = []
    consensusdir = os.path.join(dd.get_work_dir(data), "consensus")
    utils.safe_makedir(consensusdir)
    summitfiles = []
    for data in dd.sample_data_iterator(samples):
        if dd.get_chip_method(data) == "chip":
            for fn in tz.get_in(("peaks_files", "macs2"), data, []):
                if "summit" in fn:
                    summitfiles.append(fn)
                    break
        elif dd.get_chip_method(data) == "atac":
            for fn in tz.get_in(("peaks_files", "NF", "macs2") , data, []):
                if "summit" in fn:
                    summitfiles.append(fn)
    consensusfile = os.path.join(consensusdir, "consensus.bed")
    if not summitfiles:
        logger.info("No suitable peak files found, skipping consensus peak calling.")
        return samples
    consensusfile = consensus(summitfiles, consensusfile, data)
    for data in dd.sample_data_iterator(samples):
        new_samples.append([tz.assoc_in(data, ("peaks_files", "consensus"), {"main": consensusfile})])
    return new_samples
