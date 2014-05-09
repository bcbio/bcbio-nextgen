"""High level entry point for processing a sample.

Samples may include multiple lanes, or barcoded subsections of lanes,
processed together.
"""
import copy
import os

from bcbio import utils
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

    Handles merging of associated split read and discordant files if present.
    """
    if data.get("combine"):
        assert len(data["combine"].keys()) == 1
        file_key = data["combine"].keys()[0]
        extras = []
        for x in data["combine"][file_key].get("extras", []):
            if isinstance(x, (list, tuple)):
                extras.extend(x)
            else:
                extras.append(x)
        if file_key in data:
            extras.append(data[file_key])
        in_files = sorted(list(set(extras)))
        out_file = data["combine"][file_key]["out"]
        sup_exts = data.get(file_key + "-plus", {}).keys()
        for ext in sup_exts + [""]:
            merged_file = None
            if os.path.exists(utils.append_stem(out_file, "-" + ext)):
                cur_out_file, cur_in_files = out_file, []
            if ext:
                cur_in_files = list(filter(os.path.exists, (utils.append_stem(f, "-" + ext) for f in in_files)))
                cur_out_file = utils.append_stem(out_file, "-" + ext) if len(cur_in_files) > 0 else None
            else:
                cur_in_files, cur_out_file = in_files, out_file
            if cur_out_file:
                config = copy.deepcopy(data["config"])
                config["algorithm"]["save_diskspace"] = False
                if len(cur_in_files) > 0:
                    merged_file = merge_bam_files(cur_in_files, os.path.dirname(cur_out_file), config,
                                                  out_file=cur_out_file)
                else:
                    assert os.path.exists(cur_out_file)
                    merged_file = cur_out_file
            if merged_file:
                if ext:
                    data[file_key + "-plus"][ext] = merged_file
                else:
                    data[file_key] = merged_file
        data.pop("region", None)
        data.pop("combine", None)
    return [[data]]
