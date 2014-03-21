"""Retrieve run information describing files to process in a pipeline.

This handles two methods of getting processing information: from a Galaxy
next gen LIMS system or an on-file YAML configuration.
"""
import copy
import itertools
import os
from contextlib import closing

import yaml

from bcbio import utils
from bcbio.log import logger
from bcbio.galaxy.api import GalaxyApiAccess
from bcbio.illumina import flowcell
from bcbio.pipeline import alignment, config_utils, genome
from bcbio.variation import effects, genotype, population
from bcbio.variation.cortex import get_sample_name
from bcbio.bam.fastq import open_fastq

def organize(dirs, config, run_info_yaml):
    """Organize run information from a passed YAML file or the Galaxy API.

    Creates the high level structure used for subsequent processing.
    """
    if run_info_yaml and os.path.exists(run_info_yaml):
        logger.info("Using input YAML configuration: %s" % run_info_yaml)
        run_details = _run_info_from_yaml(dirs["flowcell"], run_info_yaml, config)
    else:
        logger.info("Fetching run details from Galaxy instance")
        fc_name, fc_date = flowcell.parse_dirname(dirs["flowcell"])
        galaxy_api = GalaxyApiAccess(config['galaxy_url'], config['galaxy_api_key'])
        run_details = []
        galaxy_info = galaxy_api.run_details(fc_name, fc_date)
        for item in galaxy_info["details"]:
            item["upload"] = {"method": "galaxy", "run_id": galaxy_info["run_id"],
                              "fc_name": fc_name, "fc_date": fc_date}
            run_details.append(item)
    out = []
    for item in run_details:
        # add algorithm details to configuration, avoid double specification
        item["config"] = config_utils.update_w_custom(config, item)
        item.pop("algorithm", None)
        item["dirs"] = dirs
        if "name" not in item:
            item["name"] = ["", item["description"]]
        item = add_reference_resources(item)
        out.append(item)
    return out

# ## Genome reference information

def add_reference_resources(data):
    """Add genome reference information to the item to process.
    """
    aligner = data["config"]["algorithm"].get("aligner", None)
    data["reference"] = genome.get_refs(data["genome_build"], aligner, data["dirs"]["galaxy"])
    # back compatible `sam_ref` target
    data["sam_ref"] = utils.get_in(data, ("reference", "fasta", "base"))
    ref_loc = utils.get_in(data, ("config", "resources", "species", "dir"),
                           utils.get_in(data, ("reference", "fasta", "base")))
    data["genome_resources"] = genome.get_resources(data["genome_build"], ref_loc)
    data["reference"]["snpeff"] = effects.get_snpeff_files(data)
    alt_genome = utils.get_in(data, ("config", "algorithm", "validate_genome_build"))
    if alt_genome:
        data["reference"]["alt"] = {alt_genome:
                                    genome.get_refs(alt_genome, None, data["dirs"]["galaxy"])["fasta"]}
    # Re-enable when we have ability to re-define gemini configuration directory
    if False:
        if population.do_db_build([data], check_gemini=False, need_bam=False):
            data["reference"]["gemini"] = population.get_gemini_files(data)
    return data

# ## Sample and BAM read group naming

def _clean_characters(x):
    """Clean problem characters in sample lane or descriptions.
    """
    for problem in [" ", "."]:
        x = x.replace(problem, "_")
    return x

def prep_rg_names(item, config, fc_name, fc_date):
    """Generate read group names from item inputs.
    """
    if fc_name and fc_date:
        lane_name = "%s_%s_%s" % (item["lane"], fc_date, fc_name)
    else:
        lane_name = item["description"]
    return {"rg": item["lane"],
            "sample": item["description"],
            "lane": lane_name,
            "pl": item.get("algorithm", {}).get("platform",
                                                config.get("algorithm", {}).get("platform", "illumina")).lower(),
            "pu": lane_name}

# ## Configuration file validation

def _check_for_duplicates(xs, attr, check_fn=None):
    """Identify and raise errors on duplicate items.
    """
    dups = []
    for key, vals in itertools.groupby(x[attr] for x in xs):
        if len(list(vals)) > 1:
            dups.append(key)
    if len(dups) > 0:
        psamples = []
        for x in xs:
            if x[attr] in dups:
                psamples.append(x)
        # option to skip problem based on custom input function.
        if check_fn and check_fn(psamples):
            return
        descrs = [x["description"] for x in psamples]
        raise ValueError("Duplicate '%s' found in input sample configuration.\n"
                         "Required to be unique for a project: %s\n"
                         "Problem found in these samples: %s" % (attr, dups, descrs))

def _check_for_batch_clashes(xs):
    """Check that batch names do not overlap with sample names.
    """
    names = set([x["description"] for x in xs])
    dups = set([])
    for x in xs:
        batches = utils.get_in(x, ("metadata", "batch"))
        if batches:
            if not isinstance(batches, (list, tuple)):
                batches = [batches]
            for batch in batches:
                if batch in names:
                    dups.add(batch)
    if len(dups) > 0:
        raise ValueError("Batch names must be unique from sample descriptions.\n"
                         "Clashing batch names: %s" % sorted(list(dups)))

def _check_for_misplaced(xs, subkey, other_keys):
    """Ensure configuration keys are not incorrectly nested under other keys.
    """
    problems = []
    for x in xs:
        check_dict = x.get(subkey, {})
        for to_check in other_keys:
            if to_check in check_dict:
                problems.append((x["description"], to_check, subkey))
    if len(problems) > 0:
        raise ValueError("\n".join(["Incorrectly nested keys found in sample YAML. These should be top level:",
                                    " sample         |   key name      |   nested under ",
                                    "----------------+-----------------+----------------"] +
                                   ["% 15s | % 15s | % 15s" % (a, b, c) for (a, b, c) in problems]))

ALGORITHM_KEYS = set(["platform", "aligner", "bam_clean", "bam_sort",
                      "trim_reads", "adapters", "custom_trim",
                      "align_split_size", "quality_bin",
                      "quality_format", "write_summary",
                      "merge_bamprep", "coverage", "coverage_bigwig",
                      "coverage_interval", "ploidy",
                      "variantcaller", "variant_regions",
                      "mark_duplicates", "svcaller", "recalibrate",
                      "realign", "phasing", "validate",
                      "validate_regions", "validate_genome_build",
                      "clinical_reporting", "nomap_split_size",
                      "nomap_split_targets", "ensemble", "background",
                      "disambiguate", "strandedness", "fusion_mode", "min_read_length",
                      "coverage_depth_min", "coverage_depth_max"] +
                     # back compatibility
                      ["coverage_depth"])

def _check_algorithm_keys(item):
    """Check for unexpected keys in the algorithm section.

    Needs to be manually updated when introducing new keys, but avoids silent bugs
    with typos in key names.
    """
    url = "https://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#algorithm-parameters"
    problem_keys = [k for k in item["algorithm"].iterkeys() if k not in ALGORITHM_KEYS]
    if len(problem_keys) > 0:
        raise ValueError("Unexpected configuration keyword in 'algorithm' section: %s\n"
                         "See configuration documentation for supported options:\n%s\n"
                         % (problem_keys, url))


def _detect_fastq_format(in_file, MAX_RECORDS=1000):
    ranges = {"sanger": (33, 73),
              "solexa": (59, 104),
              "illumina_1.3+": (64, 104),
              "illumina_1.5+": (66, 104),
              "illumina_1.8+": (35, 74)}

    gmin, gmax = 99, 0
    possible = set(ranges.keys())

    with closing(open_fastq(in_file)) as in_handle:
        four = itertools.islice(in_handle, 3, None, 4)
        count = 0
        for line in four:
            if len(possible) == 1:
                return possible
            if count > MAX_RECORDS:
                break
            count += 1
            vals = [ord(c) for c in line.rstrip()]
            lmin = min(vals)
            lmax = max(vals)
            for encoding, (emin, emax) in ranges.items():
                if encoding in possible:
                    if lmin < emin or lmax > emax:
                        possible.remove(encoding)

    return possible

def _check_quality_format(items):
    """
    Check if quality_format="standard" and fastq_format is not sanger
    """
    SAMPLE_FORMAT = {"illumina_1.3+": "illumina",
                     "illumina_1.5+": "illumina",
                     "illumina_1.8+": "standard",
                     "solexa": "solexa",
                     "sanger": "standard"}
    fastq_extensions = ["fq.gz", "fastq.gz", ".fastq" ".fq"]

    for item in items:
        specified_format = item["algorithm"].get("quality_format", "").lower()
        fastq_file = next((file for file in item.get('files', []) if
                           any([ext for ext in fastq_extensions if ext in file])), None)

        if fastq_file and specified_format:
            fastq_format = _detect_fastq_format(fastq_file)
            detected_encodings = set([SAMPLE_FORMAT[x] for x in fastq_format])
            if detected_encodings:
                if specified_format not in detected_encodings:
                    raise ValueError("Quality format specified in the YAML "
                                     "file might be a different encoding. "
                                     "'%s' was specified but possible formats "
                                     "detected were %s." % (specified_format,
                                                            ", ".join(detected_encodings)))


def _check_aligner(item):
    """Ensure specified aligner is valid choice.
    """
    allowed = set(alignment.TOOLS.keys() + [None, False])
    if item["algorithm"].get("aligner") not in allowed:
        raise ValueError("Unexpected algorithm 'aligner' parameter: %s\n"
                         "Supported options: %s\n" %
                         (item["algorithm"].get("aligner"), sorted(list(allowed))))

def _check_variantcaller(item):
    """Ensure specified variantcaller is a valid choice.
    """
    allowed = set(genotype.get_variantcallers().keys() + [None, False])
    vcs = item["algorithm"].get("variantcaller", "gatk")
    if not isinstance(vcs, (tuple, list)):
        vcs = [vcs]
    problem = [x for x in vcs if x not in allowed]
    if len(problem) > 0:
        raise ValueError("Unexpected algorithm 'variantcaller' parameter: %s\n"
                         "Supported options: %s\n" % (problem, sorted(list(allowed))))

def _check_sample_config(items, in_file):
    """Identify common problems in input sample configuration files.
    """
    logger.info("Checking sample YAML configuration: %s" % in_file)
    _check_quality_format(items)
    _check_for_duplicates(items, "lane")
    _check_for_duplicates(items, "description")
    _check_for_batch_clashes(items)
    _check_for_misplaced(items, "algorithm",
                         ["resources", "metadata", "analysis",
                          "description", "genome_build", "lane", "files"])

    [_check_algorithm_keys(x) for x in items]
    [_check_aligner(x) for x in items]
    [_check_variantcaller(x) for x in items]


# ## Read bcbio_sample.yaml files

def _normalize_files(item, fc_dir):
    """Ensure the files argument is a list of absolute file names.
    Handles BAM, single and paired end fastq.
    """
    files = item.get("files")
    if files:
        if isinstance(files, basestring):
            files = [files]
        if fc_dir:
            fastq_dir = flowcell.get_fastq_dir(fc_dir)
        else:
            fastq_dir = os.getcwd()
        files = [x if os.path.isabs(x) else os.path.normpath(os.path.join(fastq_dir, x))
                 for x in files]
        _sanity_check_files(item, files)
        item["files"] = files
    return item

def _sanity_check_files(item, files):
    """Ensure input files correspond with supported
    """
    msg = None
    file_types = set([("bam" if x.endswith(".bam") else "fastq") for x in files if x])
    if len(file_types) > 1:
        msg = "Found multiple file types (BAM and fastq)"
    file_type = file_types.pop()
    if file_type == "bam":
        if len(files) != 1:
            msg = "Expect a single BAM file input as input"
    elif file_type == "fastq":
        if len(files) not in [1, 2]:
            msg = "Expect either 1 (single end) or 2 (paired end) fastq inputs"
    if msg:
        raise ValueError("%s for %s: %s" % (msg, item.get("description", ""), files))

def _run_info_from_yaml(fc_dir, run_info_yaml, config):
    """Read run information from a passed YAML file.
    """
    with open(run_info_yaml) as in_handle:
        loaded = yaml.load(in_handle)
    fc_name, fc_date = None, None
    if fc_dir:
        try:
            fc_name, fc_date = flowcell.parse_dirname(fc_dir)
        except ValueError:
            pass
    global_config = {}
    global_vars = {}
    if isinstance(loaded, dict):
        global_config = copy.deepcopy(loaded)
        del global_config["details"]
        if "fc_name" in loaded and "fc_date" in loaded:
            fc_name = loaded["fc_name"].replace(" ", "_")
            fc_date = str(loaded["fc_date"]).replace(" ", "_")
        global_vars = global_config.pop("globals", {})
        loaded = loaded["details"]

    run_details = []
    for i, item in enumerate(loaded):
        item = _normalize_files(item, fc_dir)
        if "lane" not in item:
            item["lane"] = str(i + 1)
        item["lane"] = _clean_characters(str(item["lane"]))
        if "description" not in item:
            if len(item.get("files", [])) == 1 and item["files"][0].endswith(".bam"):
                item["description"] = get_sample_name(item["files"][0])
            else:
                raise ValueError("No `description` sample name provided for input #%s" % (i + 1))
        item["description"] = _clean_characters(str(item["description"]))
        upload = global_config.get("upload", {})
        # Handle specifying a local directory directly in upload
        if isinstance(upload, basestring):
            upload = {"dir": upload}
        if fc_name and fc_date:
            upload["fc_name"] = fc_name
            upload["fc_date"] = fc_date
        upload["run_id"] = ""
        item["upload"] = upload
        item["algorithm"] = _replace_global_vars(item["algorithm"], global_vars)
        item["algorithm"] = genome.abs_file_paths(item["algorithm"],
                                                  ignore_keys=["variantcaller", "realign", "recalibrate",
                                                               "phasing", "svcaller"])

        item["rgnames"] = prep_rg_names(item, config, fc_name, fc_date)
        item["test_run"] = global_config.get("test_run", False)
        run_details.append(item)
    _check_sample_config(run_details, run_info_yaml)
    return run_details

def _replace_global_vars(xs, global_vars):
    """Replace globally shared names from input header with value.

    The value of the `algorithm` item may be a pointer to a real
    file specified in the `global` section. If found, replace with
    the full value.
    """
    if isinstance(xs, (list, tuple)):
        return [_replace_global_vars(x) for x in xs]
    elif isinstance(xs, dict):
        final = {}
        for k, v in xs.iteritems():
            if isinstance(v, basestring) and v in global_vars:
                v = global_vars[v]
            final[k] = v
        return final
    else:
        return xs
