"""High level entry point for processing a sample.

Samples may include multiple lanes, or barcoded subsections of lanes,
processed together.
"""
import os
import subprocess


from bcbio.utils import file_exists, save_diskspace
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline.merge import (combine_fastq_files, merge_bam_files)
from bcbio.pipeline.qcsummary import generate_align_summary
from bcbio.rnaseq.cufflinks import assemble_transcripts
from bcbio.pipeline import shared
from bcbio.rnaseq import count

def merge_sample(data):
    """Merge fastq and BAM files for multiple samples.
    """
    logger.info("Combining fastq and BAM files %s" % str(data["name"]))
    config = shared.update_config_w_custom(data["config"], data["info"])
    genome_build, sam_ref = shared.ref_genome_info(data["info"], config, data["dirs"])
    if config["algorithm"].get("upload_fastq", False):
        fastq1, fastq2 = combine_fastq_files(data["fastq_files"], data["dirs"]["work"],
                                             config)
    else:
        fastq1, fastq2 = None, None
    sort_bam = merge_bam_files(data["bam_files"], data["dirs"]["work"], config)
    return [[{"name": data["name"], "metadata": data["info"].get("metadata", {}),
              "info": data["info"],
              "genome_build": genome_build, "sam_ref": sam_ref,
              "work_bam": sort_bam, "fastq1": fastq1, "fastq2": fastq2,
              "dirs": data["dirs"], "config": config,
              "config_file": data["config_file"]}]]

# ## General processing

def process_sample(data):
    """Finalize processing for a sample, potentially multiplexed.
    """
    if data["config"]["algorithm"].get("transcript_assemble", False):
        data["tx_file"] = assemble_transcripts(data["work_bam"], data["sam_ref"],
                                               data["config"])
    if data["sam_ref"] is not None:
        logger.info("Generating summary files: %s" % str(data["name"]))
        data["summary"] = generate_align_summary(data["work_bam"], data["fastq2"] is not None,
                                                 data["sam_ref"], data["name"],
                                                 data["config"], data["dirs"])
    return [[data]]

def generate_transcript_counts(data):
    """Generate counts per transcript from an alignment"""
    data["count_file"] = count.htseq_count(data)
    return [[data]]

def generate_bigwig(data):
    """Provide a BigWig coverage file of the sorted alignments.
    """
    logger.info("Preparing BigWig file %s" % str(data["name"]))
    bam_file = data["work_bam"]
    wig_file = "%s.bigwig" % os.path.splitext(bam_file)[0]
    if not file_exists(wig_file):
        with file_transaction(wig_file) as tx_file:
            cl = ["bam_to_wiggle.py", bam_file,
                  data["config_file"], "--outfile=%s" % tx_file]
            subprocess.check_call(cl)
    return [[data]]
