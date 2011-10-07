"""High level entry point for processing a sample.

Samples may include multiple lanes, or barcoded subsections of lanes,
processed together.
"""
import os
import subprocess


from bcbio.utils import file_transaction
from bcbio.pipeline.lane import _update_config_w_custom
from bcbio.pipeline import log
from bcbio.pipeline.merge import (combine_fastq_files, merge_bam_files)
from bcbio.pipeline.qcsummary import generate_align_summary
from bcbio.pipeline.variation import (recalibrate_quality, finalize_genotyper,
                                      variation_effects)
from bcbio.rnaseq.cufflinks import assemble_transcripts
from bcbio.pipeline.shared import ref_genome_info

def merge_sample(data):
    """Merge fastq and BAM files for multiple samples.
    """
    log.info("Combining fastq and BAM files %s" % str(data["name"]))
    config = _update_config_w_custom(data["config"], data["info"])
    genome_build, sam_ref = ref_genome_info(data["info"], config, data["dirs"])
    fastq1, fastq2 = combine_fastq_files(data["fastq_files"], data["dirs"]["work"])
    sort_bam = merge_bam_files(data["bam_files"], data["dirs"]["work"], config)
    return [[{"name": data["name"],
              "genome_build": genome_build, "sam_ref": sam_ref,
              "work_bam": sort_bam, "fastq1": fastq1, "fastq2": fastq2,
              "dirs": data["dirs"], "config": config,
              "config_file": data["config_file"]}]]

def recalibrate_sample(data):
    """Recalibrate quality values from aligned sample BAM file.
    """
    log.info("Recalibrating %s with GATK" % str(data["name"]))
    if data["config"]["algorithm"]["recalibrate"]:
        data["work_bam"] = recalibrate_quality(data["work_bam"], data["fastq1"],
                                               data["fastq2"], data["sam_ref"],
                                               data["dirs"], data["config"])
    return [[data]]

# ## General processing

def process_sample(data):
    """Finalize processing for a sample, potentially multiplexed.
    """
    if data["config"]["algorithm"]["snpcall"]:
        log.info("Finalizing variant calls %s with GATK" % str(data["name"]))
        data["vrn_file"] = finalize_genotyper(data["vrn_file"], data["sam_ref"],
                                              data["config"])
        log.info("Calculating variation effects for %s" % str(data["name"]))
        ann_vrn_file, effects_file = variation_effects(data["vrn_file"],
                                                       data["sam_ref"],
                                                       data["genome_build"],
                                                       data["config"])
        if ann_vrn_file:
            data["vrn_file"] = ann_vrn_file
            data["effects_file"] = effects_file
    if data["config"]["algorithm"].get("transcript_assemble", False):
        data["tx_file"] = assemble_transcripts(data["work_bam"], data["sam_ref"],
                                               data["config"])
    if data["sam_ref"] is not None:
        log.info("Generating summary files: %s" % str(data["name"]))
        generate_align_summary(data["work_bam"], data["fastq2"] is not None,
                               data["sam_ref"], data["name"],
                               data["config"], data["dirs"])
    log.info("Preparing BigWig file %s" % str(data["name"]))
    data["bigwig_file"] = bam_to_wig(data["work_bam"], data["config"],
                                     data["config_file"])
    return [[data]]

def bam_to_wig(bam_file, config, config_file):
    """Provide a BigWig coverage file of the sorted alignments.
    """
    wig_file = "%s.bigwig" % os.path.splitext(bam_file)[0]
    if not (os.path.exists(wig_file) and os.path.getsize(wig_file) > 0):
        cl = [config["analysis"]["towig_script"], bam_file, config_file]
        with file_transaction(wig_file):
            subprocess.check_call(cl)
    return wig_file

