"""Examine sequencing coverage, identifying transcripts lacking sufficient coverage for variant calling.

Handles identification of low coverage regions in defined genes of interest or the entire transcript.
"""
import collections
import copy
import os

import yaml

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import config_utils
from bcbio.provenance import do

def _prep_coverage_file(species, covdir, config):
    """Ensure input coverage file is correct, handling special keywords for whole exome.
    Returns the input coverage file and keyword for special cases.
    """
    cov_file = config["algorithm"]["coverage"]
    cov_kw = None
    if cov_file == "exons":
        cov_kw = cov_file
        cov_file = os.path.join(covdir, "%s-%s.txt" % (species, cov_kw))
    else:
        cov_file = os.path.normpath(os.path.join(os.path.split(covdir)[0], cov_file))
        assert os.path.exists(cov_file), \
            "Did not find input file for coverage: %s" % cov_file
    return cov_file, cov_kw

def _prep_coverage_config(samples, config):
    """Create input YAML configuration and directories for running coverage assessment.
    """
    covdir = utils.safe_makedir(os.path.join(samples[0]["dirs"]["work"], "coverage"))
    name = samples[0]["name"][-1].replace(" ", "_")
    cur_covdir = utils.safe_makedir(os.path.join(covdir, name))
    out_file = os.path.join(cur_covdir, "coverage_summary.csv")
    config_file = os.path.join(cur_covdir, "coverage-in.yaml")
    species = samples[0]["genome_resources"]["aliases"]["ensembl"]
    cov_file, cov_kw = _prep_coverage_file(species, covdir, config)
    out = {"params": {"species": species,
                      "build": samples[0]["genome_build"],
                      "coverage": 13,
                      "transcripts": "canonical",
                      "block": {"min": 100, "distance": 10}},
           "regions": cov_file,
           "ref-file": os.path.abspath(samples[0]["sam_ref"]),
           "experiments": []
           }
    if cov_kw:
        out["params"]["regions"] = cov_kw
    for data in samples:
        out["experiments"].append({"name": data["name"][-1],
                                   "samples": [{"coverage": str(data["work_bam"])}]})
    with open(config_file, "w") as out_handle:
        yaml.dump(out, out_handle, allow_unicode=False, default_flow_style=False)
    return config_file, out_file

def summary(samples, config):
    """Provide summary information on a single sample across regions of interest.
    """
    try:
        bc_jar = config_utils.get_jar("bcbio.coverage", config_utils.get_program("bcbio_coverage", config, "dir"))
    except ValueError:
        logger.warning("No coverage calculations: Did not find bcbio.coverage jar from system config")
        return [[x] for x in samples]
    config_file, out_file = _prep_coverage_config(samples, config)
    tmp_dir = utils.safe_makedir(os.path.join(os.path.dirname(out_file), "tmp"))
    resources = config_utils.get_resources("bcbio_coverage", config)
    config = copy.deepcopy(config)
    config["algorithm"]["memory_adjust"] = {"direction": "increase",
                                            "magnitude": config["algorithm"].get("num_cores", 1)}
    jvm_opts = config_utils.adjust_opts(resources.get("jvm_opts", ["-Xms750m", "-Xmx2g"]), config)
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            java_args = ["-Djava.io.tmpdir=%s" % tmp_dir, "-Djava.awt.headless=true"]
            cmd = ["java"] + jvm_opts + java_args + ["-jar", bc_jar, "multicompare", config_file,
                                                     tx_out_file, "-c", str(config["algorithm"].get("num_cores", 1))]
            do.run(cmd, "Summarizing coverage with bcbio.coverage", samples[0])
    out = []
    for x in samples:
        x["coverage"] = {"summary": out_file}
        out.append([x])
    return out

def summarize_samples(samples, run_parallel):
    """Provide summary information for sample coverage across regions of interest.
    """
    to_run = collections.defaultdict(list)
    extras = []
    for data in [x[0] for x in samples]:
        if ("coverage" in data["config"]["algorithm"] and
              data["genome_resources"].get("aliases", {}).get("ensembl")):
            to_run[(data["genome_build"], data["config"]["algorithm"]["coverage"])].append(data)
        else:
            extras.append([data])
    out = []
    if len(to_run) > 0:
        args = []
        for sample_group in to_run.itervalues():
            config = sample_group[0]["config"]
            args.append((sample_group, config))
        out.extend(run_parallel("coverage_summary", args))
    return out + extras
