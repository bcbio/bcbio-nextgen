"""Integration with Galaxy nglims.
"""
import collections
import copy
import glob
import operator
import os
import string
import subprocess

import joblib
import yaml

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.galaxy.api import GalaxyApiAccess
from bcbio.illumina import flowcell

def prep_samples_and_config(run_folder, ldetails, fastq_dir, config):
    """Prepare sample fastq files and provide global sample configuration for the flowcell.

    Handles merging of fastq files split by lane and also by the bcl2fastq
    preparation process.
    """
    fastq_final_dir = utils.safe_makedir(os.path.join(fastq_dir, "merged"))
    cores = utils.get_in(config, ("algorithm", "num_cores"), 1)
    ldetails = joblib.Parallel(cores)(joblib.delayed(_prep_sample_and_config)(x, fastq_dir, fastq_final_dir)
                                      for x in _group_same_samples(ldetails))
    config_file = _write_sample_config(run_folder, ldetails)
    return config_file, fastq_final_dir

def _prep_sample_and_config(ldetail_group, fastq_dir, fastq_final_dir):
    """Prepare output fastq file and configuration for a single sample.
    """
    files = []
    print "->", ldetail_group[0]["name"], len(ldetail_group)
    for read in ["R1", "R2"]:
        fastq_inputs = reduce(operator.add, (_get_fastq_files(x, read, fastq_dir) for x in ldetail_group))
        if len(fastq_inputs) > 0:
            files.append(_concat_bgzip_fastq(fastq_inputs, fastq_final_dir, read, ldetail_group[0]))
    if len(files) == 0:
        raise ValueError("Did not find input fastq files for %s %s" % (ldetail["project_name"], ldetail["name"]))
    out = ldetail_group[0]
    out["files"] = files
    return out

def _write_sample_config(run_folder, ldetails):
    """Generate a bcbio-nextgen YAML configuration file for processing a sample.

    TODO: Fill in reasonable defaults for algorithm configuration from default templates.
    """
    out_file = os.path.join(run_folder, "%s.yaml" % os.path.basename(run_folder))
    with open(out_file, "w") as out_handle:
        fc_name, fc_date = flowcell.parse_dirname(run_folder)
        out = {"details": sorted([_prepare_sample(x, run_folder) for x in ldetails],
                                 key=operator.itemgetter("name", "description")),
               "fc_name": fc_name,
               "fc_date": fc_date}
        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out_file

def _prepare_sample(data, run_folder):
    """Extract passed keywords from input LIMS information.
    """
    want = set(["description", "files", "genome_build", "name", "analysis", "upload", "algorithm"])
    out = {}
    for k, v in data.items():
        if k in want:
            out[k] = _relative_paths(v, run_folder)
    return out

def _relative_paths(xs, base_path):
    """Adjust paths to be relative to the provided base path.
    """
    if isinstance(xs, basestring):
        if xs.startswith(base_path):
            return xs.replace(base_path + "/", "", 1)
        else:
            return xs
    elif isinstance(xs, (list, tuple)):
        return [_relative_paths(x, base_path) for x in xs]
    elif isinstance(xs, dict):
        out = {}
        for k, v in xs.items():
            out[k] = _relative_paths(v, base_path)
        return out
    else:
        return xs

def _get_fastq_files(ldetail, read, fastq_dir):
    """Retrieve fastq files corresponding to the sample and read number.
    """
    return glob.glob(os.path.join(fastq_dir, "Project_%s" % ldetail["project_name"],
                                  "Sample_%s" % ldetail["name"],
                                  "%s_*_%s_*.fastq.gz" % (ldetail["name"], read)))

def _concat_bgzip_fastq(finputs, out_dir, read, ldetail):
    """Concatenate multiple input fastq files, preparing a bgzipped output file.
    """
    out_file = os.path.join(out_dir, "%s_%s.fastq.gz" % (ldetail["name"], read))
    if False:
    # if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            subprocess.check_call("zcat %s | bgzip -c > %s" % (" ".join(finputs), tx_out_file), shell=True)
    return out_file

def _group_same_samples(ldetails):
    """Move samples into groups -- same groups have identical projects and descriptions
    """
    sample_groups = collections.defaultdict(list)
    for ldetail in ldetails:
        sample_groups[(ldetail["project_name"], ldetail["description"])].append(ldetail)
    return sample_groups.values()

def get_runinfo(galaxy_url, galaxy_apikey, run_folder):
    """Retrieve flattened run information for a processed directory from Galaxy nglims API.
    """
    galaxy_api = GalaxyApiAccess(galaxy_url, galaxy_apikey)
    fc_name, fc_date = flowcell.parse_dirname(run_folder)
    galaxy_info = galaxy_api.run_details(fc_name, fc_date)
    ldetails = _flatten_lane_details(galaxy_info)
    out = []
    for item in ldetails:
        item["upload"] = {"method": "galaxy", "run_id": galaxy_info["run_id"],
                          "fc_name": fc_name, "fc_date": fc_date,
                          "galaxy_url": galaxy_url, "galaxy_api_key": galaxy_apikey}
        for k in ["lab_association", "private_libs", "researcher", "researcher_id", "sample_id",
                  "galaxy_library", "galaxy_role"]:
            item["upload"][k] = item.pop(k, "")
        out.append(item)
    return out

def _flatten_lane_details(runinfo):
    """Provide flattened lane information with multiplexed barcodes separated.
    """
    out = []
    for ldetail in runinfo["details"]:
        for i, barcode in enumerate(ldetail.get("multiplex", [{}])):
            cur = copy.deepcopy(ldetail)
            cur["name"] = "%s-%s" % (ldetail["name"], i + 1)
            cur["description"] = barcode.get("name", ldetail["description"])
            cur["bc_index"] = barcode.get("sequence", "")
            cur["project_name"] = clean_name(ldetail["project_name"])
            out.append(cur)
    return out

def clean_name(xs):
    final = []
    safec = "_"
    for x in xs:
        if x not in string.ascii_letters + string.digits:
            if len(final) > 0 and final[-1] != safec:
                final.append(safec)
        else:
            final.append(x)
    if final[-1] == safec:
        final = final[:-1]
    return "".join(final)
