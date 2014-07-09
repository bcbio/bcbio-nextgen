"""High level entry point for processing a sample.

Samples may include multiple lanes, or barcoded subsections of lanes,
processed together.
"""
import os
import copy
import subprocess

from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline.merge import (combine_fastq_files, merge_bam_files)
from bcbio.pipeline import config_utils

# ## Merging

def merge_sample(data):
    """Merge fastq and BAM files for multiple samples.
    """
    logger.debug("Combining fastq and BAM files %s" % str(data["name"]))
    config = config_utils.update_w_custom(data["config"], data["info"])
    if config["algorithm"].get("upload_fastq", False):
        fastq1, fastq2 = combine_fastq_files(data["fastq_files"], data["dirs"]["work"],
                                             config)
    else:
        fastq1, fastq2 = None, None

    out_file = os.path.join(data["dirs"]["work"],
                            data["info"]["rgnames"]["sample"] + ".bam")
    sort_bam = merge_bam_files(data["bam_files"], data["dirs"]["work"],
                               config, out_file=out_file)
    return [[{"name": data["name"], "metadata": data["info"].get("metadata", {}),
              "info": data["info"],
              "genome_build": data["genome_build"], "sam_ref": data["sam_ref"],
              "work_bam": sort_bam, "fastq1": fastq1, "fastq2": fastq2,
              "dirs": data["dirs"], "config": config,
              "config_file": data["config_file"]}]]

def delayed_bam_merge(data):
    """Perform a merge on previously prepped files, delayed in processing.
    """
    if data.get("combine"):
        assert len(data["combine"].keys()) == 1
        file_key = data["combine"].keys()[0]
        in_files = sorted(list(set([data[file_key]] + data["combine"][file_key].get("extras", []))))
        out_file = data["combine"][file_key]["out"]
        logger.debug("Combining BAM files to %s" % out_file)
        config = copy.deepcopy(data["config"])
        config["algorithm"]["save_diskspace"] = False
        merged_file = merge_bam_files(in_files, os.path.dirname(out_file), config,
                                      out_file=out_file)
        if data.has_key("region"):
            del data["region"]
        del data["combine"]
        data[file_key] = merged_file
    return [[data]]

# ## General processing

def generate_bigwig(data):
    """Provide a BigWig coverage file of the sorted alignments.
    """
    if data["config"]["algorithm"].get("coverage_bigwig", True):
        logger.info("Preparing BigWig file %s" % str(data["name"]))
        bam_file = data["work_bam"]
        wig_file = "%s.bigwig" % os.path.splitext(bam_file)[0]
        if not file_exists(wig_file):
            with file_transaction(wig_file) as tx_file:
                cl = ["bam_to_wiggle.py", bam_file,
                      data["config_file"], "--outfile=%s" % tx_file]
                subprocess.check_call(cl)
    return [[data]]
