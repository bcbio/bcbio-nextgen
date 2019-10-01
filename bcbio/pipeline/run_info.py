"""Retrieve run information describing files to process in a pipeline.

This handles two methods of getting processing information: from a Galaxy
next gen LIMS system or an on-file YAML configuration.
"""
import collections
from contextlib import closing
import copy
import glob
import itertools
import operator
import os
import string

import six
import toolz as tz
import yaml
from bcbio import install, utils, structural
from bcbio.bam import fastq, ref
from bcbio.log import logger
from bcbio.distributed import objectstore
from bcbio.illumina import flowcell
from bcbio.pipeline import alignment, config_utils, genome
from bcbio.pipeline import datadict as dd
from bcbio.provenance import diagnostics, programs, versioncheck
from bcbio.provenance import data as provenancedata
from bcbio.qc import viral

from bcbio.variation import annotation, effects, genotype, population, joint, vcfutils, vcfanno
from bcbio.variation.cortex import get_sample_name
from bcbio.bam.fastq import open_fastq
from functools import reduce

ALLOWED_CONTIG_NAME_CHARS = set(list(string.digits) + list(string.ascii_letters) + ["-", "_", "*", ":", "."])
ALGORITHM_NOPATH_KEYS = ["variantcaller", "realign", "recalibrate", "peakcaller",
                         "expression_caller", "singlecell_quantifier",
                         "fusion_caller",
                         "svcaller", "hetcaller", "jointcaller", "tools_off", "tools_on",
                         "mixup_check", "qc", "transcript_assembler"]
ALGORITHM_FILEONLY_KEYS = ["custom_trim", "vcfanno"]
# these analysis pipelines use R heavily downstream, and need to have samplenames
# cleaned up to conform to R specifications
R_DOWNSTREAM_ANALYSIS = ["rna-seq", "fastrna-seq", "scrna-seq", "chip-seq",
                         "scrna-seq"]

def organize(dirs, config, run_info_yaml, sample_names=None, is_cwl=False,
             integrations=None):
    """Organize run information from a passed YAML file or the Galaxy API.

    Creates the high level structure used for subsequent processing.

    sample_names is a list of samples to include from the overall file, for cases
    where we are running multiple pipelines from the same configuration file.
    """
    from bcbio.pipeline import qcsummary
    if integrations is None: integrations = {}
    logger.info("Using input YAML configuration: %s" % run_info_yaml)
    assert run_info_yaml and os.path.exists(run_info_yaml), \
        "Did not find input sample YAML file: %s" % run_info_yaml
    run_details = _run_info_from_yaml(dirs, run_info_yaml, config, sample_names,
                                      is_cwl=is_cwl, integrations=integrations)
    remote_retriever = None
    for iname, retriever in integrations.items():
        if iname in config:
            run_details = retriever.add_remotes(run_details, config[iname])
            remote_retriever = retriever
    out = []
    for item in run_details:
        item["dirs"] = dirs
        if "name" not in item:
            item["name"] = ["", item["description"]]
        elif isinstance(item["name"], six.string_types):
            description = "%s-%s" % (item["name"], clean_name(item["description"]))
            item["name"] = [item["name"], description]
            item["description"] = description
        # add algorithm details to configuration, avoid double specification
        item["resources"] = _add_remote_resources(item["resources"])
        item["config"] = config_utils.update_w_custom(config, item)
        item.pop("algorithm", None)
        item = add_reference_resources(item, remote_retriever)
        item["config"]["algorithm"]["qc"] = qcsummary.get_qc_tools(item)
        item["config"]["algorithm"]["vcfanno"] = vcfanno.find_annotations(item, remote_retriever)
        # Create temporary directories and make absolute, expanding environmental variables
        tmp_dir = tz.get_in(["config", "resources", "tmp", "dir"], item)
        if tmp_dir:
            # if no environmental variables, make and normalize the directory
            # otherwise we normalize later in distributed.transaction:
            if os.path.expandvars(tmp_dir) == tmp_dir:
                tmp_dir = utils.safe_makedir(os.path.expandvars(tmp_dir))
                tmp_dir = genome.abs_file_paths(tmp_dir, do_download=not integrations)
            item["config"]["resources"]["tmp"]["dir"] = tmp_dir
        out.append(item)
    out = _add_provenance(out, dirs, config, not is_cwl)
    return out

def normalize_world(data):
    """Normalize a data object, useful after serializetion via CWL.
    """
    data = _normalize_files(data)
    return data

def _add_provenance(items, dirs, config, add_provenance=True):
    if add_provenance:
        p = programs.write_versions(dirs, config=config)
        d = provenancedata.write_versions(dirs, items)
        versioncheck.testall(items)
        p_db = diagnostics.initialize(dirs)
    out = []
    for item in items:
        if add_provenance:
            entity_id = diagnostics.store_entity(item)
            item["config"]["resources"]["program_versions"] = p
            item["provenance"] = {"programs": p, "entity": entity_id,
                                  "db": p_db, "data": d}
        out.append([item])
    return out

def setup_directories(work_dir, fc_dir, config, config_file):
    fastq_dir, galaxy_dir, config_dir = _get_full_paths(flowcell.get_fastq_dir(fc_dir)
                                                        if fc_dir else None,
                                                        config, config_file)
    # check default install for tool data if not found locally
    if not os.path.exists(os.path.join(galaxy_dir, "tool-data")):
        _, config_file = config_utils.load_system_config(work_dir=work_dir, allow_missing=True)
        if config_file and os.path.exists(os.path.join(os.path.dirname(config_file), "tool-data")):
            galaxy_dir = os.path.dirname(config_file)
    return {"fastq": fastq_dir, "galaxy": galaxy_dir,
            "work": work_dir, "flowcell": fc_dir, "config": config_dir}

def _get_full_paths(fastq_dir, config, config_file):
    """Retrieve full paths for directories in the case of relative locations.
    """
    if fastq_dir:
        fastq_dir = utils.add_full_path(fastq_dir)
    config_dir = utils.add_full_path(os.path.dirname(config_file))
    galaxy_config_file = utils.add_full_path(config.get("galaxy_config", "universe_wsgi.ini"),
                                             config_dir)
    return fastq_dir, os.path.dirname(galaxy_config_file), config_dir

# ## Remote resources

def _add_remote_resources(resources):
    """Retrieve remote resources like GATK/MuTect jars present in S3.
    """
    out = copy.deepcopy(resources)
    for prog, info in resources.items():
        for key, val in info.items():
            if key == "jar" and objectstore.is_remote(val):
                store_dir = utils.safe_makedir(os.path.join(os.getcwd(), "inputs", "jars", prog))
                fname = objectstore.download(val, store_dir, store_dir)
                version_file = os.path.join(store_dir, "version.txt")
                if not utils.file_exists(version_file):
                    version = install.get_gatk_jar_version(prog, fname)
                    with open(version_file, "w") as out_handle:
                        out_handle.write(version)
                else:
                    with open(version_file) as in_handle:
                        version = in_handle.read().strip()
                del out[prog][key]
                out[prog]["dir"] = store_dir
                out[prog]["version"] = version
    return out

# ## Genome reference information

def add_reference_resources(data, remote_retriever=None):
    """Add genome reference information to the item to process.
    """
    aligner = data["config"]["algorithm"].get("aligner", None)
    if remote_retriever:
        data["reference"] = remote_retriever.get_refs(data["genome_build"],
                                                      alignment.get_aligner_with_aliases(aligner, data),
                                                      data["config"])
    else:
        data["reference"] = genome.get_refs(data["genome_build"], alignment.get_aligner_with_aliases(aligner, data),
                                            data["dirs"]["galaxy"], data)
        _check_ref_files(data["reference"], data)
    # back compatible `sam_ref` target
    data["sam_ref"] = utils.get_in(data, ("reference", "fasta", "base"))
    ref_loc = utils.get_in(data, ("config", "resources", "species", "dir"),
                           utils.get_in(data, ("reference", "fasta", "base")))
    if remote_retriever:
        data = remote_retriever.get_resources(data["genome_build"], ref_loc, data)
    else:
        data["genome_resources"] = genome.get_resources(data["genome_build"], ref_loc, data)
    data["genome_resources"] = genome.add_required_resources(data["genome_resources"])
    if effects.get_type(data) == "snpeff" and "snpeff" not in data["reference"]:
        data["reference"]["snpeff"] = effects.get_snpeff_files(data)
    if "genome_context" not in data["reference"]:
        data["reference"]["genome_context"] = annotation.get_context_files(data)
    if "viral" not in data["reference"]:
        data["reference"]["viral"] = viral.get_files(data)
    if not data["reference"]["viral"]:
        data["reference"]["viral"] = None
    if "versions" not in data["reference"]:
        data["reference"]["versions"] = _get_data_versions(data)

    data = _fill_validation_targets(data)
    data = _fill_prioritization_targets(data)
    data = _fill_capture_regions(data)
    # Re-enable when we have ability to re-define gemini configuration directory
    if False:
        data["reference"]["gemini"] = population.get_gemini_files(data)
    return data

def _get_data_versions(data):
    """Retrieve CSV file with version information for reference data.
    """
    genome_dir = install.get_genome_dir(data["genome_build"], data["dirs"].get("galaxy"), data)
    if genome_dir:
        version_file = os.path.join(genome_dir, "versions.csv")
        if version_file and os.path.exists(version_file):
            return version_file
    return None

def _check_ref_files(ref_info, data):
    problems = []
    if not data["genome_build"]:
        problems.append("Did not find 'genome_build' for sample: %s" % dd.get_sample_name(data))
    elif not tz.get_in(["fasta", "base"], ref_info):
        problems.append("Did not find fasta reference file for genome %s.\n" % (data["genome_build"]) +
                        "Check tool-data/*.loc files to ensure paths to reference data are correct.")
    else:
        for contig in ref.file_contigs(ref_info["fasta"]["base"], data["config"]):
            cur_problems = set([])
            for char in list(contig.name):
                if char not in ALLOWED_CONTIG_NAME_CHARS:
                    cur_problems.add(char)
            if len(cur_problems) > 0:
                problems.append("Found non-allowed characters in chromosome name %s: %s" %
                                (contig.name, " ".join(list(cur_problems))))
    if len(problems) > 0:
        msg = ("\nProblems with input reference file %s\n" % tz.get_in(["fasta", "base"], ref_info))
        raise ValueError(msg + "\n".join(problems) + "\n")

def _fill_validation_targets(data):
    """Fill validation targets pointing to globally installed truth sets.
    """
    ref_file = dd.get_ref_file(data)
    sv_truth = tz.get_in(["config", "algorithm", "svvalidate"], data, {})
    sv_targets = (zip(itertools.repeat("svvalidate"), sv_truth.keys()) if isinstance(sv_truth, dict)
                  else [["svvalidate"]])
    for vtarget in [list(xs) for xs in [["validate"], ["validate_regions"], ["variant_regions"]] + list(sv_targets)]:
        val = tz.get_in(["config", "algorithm"] + vtarget, data)
        if val and not os.path.exists(val) and not objectstore.is_remote(val):
            installed_val = os.path.normpath(os.path.join(os.path.dirname(ref_file), os.pardir, "validation", val))
            if os.path.exists(installed_val):
                data = tz.update_in(data, ["config", "algorithm"] + vtarget, lambda x: installed_val)
            else:
                raise ValueError("Configuration problem. Validation file not found for %s: %s" %
                                 (vtarget, val))
    return data

def _fill_capture_regions(data):
    """Fill short-hand specification of BED capture regions.
    """
    special_targets = {"sv_regions": ("exons", "transcripts")}
    ref_file = dd.get_ref_file(data)
    for target in ["variant_regions", "sv_regions", "coverage"]:
        val = tz.get_in(["config", "algorithm", target], data)
        if val and not os.path.exists(val) and not objectstore.is_remote(val):
            installed_vals = []
            # Check prioritize directory
            for ext in [".bed", ".bed.gz"]:
                installed_vals += glob.glob(os.path.normpath(os.path.join(os.path.dirname(ref_file), os.pardir,
                                                                          "coverage", val + ext)))
            if len(installed_vals) == 0:
                if target not in special_targets or not val.startswith(special_targets[target]):
                    raise ValueError("Configuration problem. BED file not found for %s: %s" %
                                     (target, val))
            else:
                assert len(installed_vals) == 1, installed_vals
                data = tz.update_in(data, ["config", "algorithm", target], lambda x: installed_vals[0])
    return data

def _fill_prioritization_targets(data):
    """Fill in globally installed files for prioritization.
    """
    ref_file = dd.get_ref_file(data)
    for target in ["svprioritize", "coverage"]:
        val = tz.get_in(["config", "algorithm", target], data)
        if val and not os.path.exists(val) and not objectstore.is_remote(val):
            installed_vals = []
            # Check prioritize directory
            for ext in [".bed", ".bed.gz"]:
                installed_vals += glob.glob(os.path.normpath(os.path.join(os.path.dirname(ref_file), os.pardir,
                                                                          "coverage", "prioritize",
                                                                          val + "*%s" % ext)))
            # Check sv-annotation directory for prioritize gene name lists
            if target == "svprioritize":
                simple_sv_bin = utils.which("simple_sv_annotation.py")
                if simple_sv_bin:
                    installed_vals += glob.glob(os.path.join(os.path.dirname(os.path.realpath(simple_sv_bin)),
                                                             "%s*" % os.path.basename(val)))
            if len(installed_vals) == 0:
                # some targets can be filled in later
                if target not in set(["coverage"]):
                    raise ValueError("Configuration problem. BED file not found for %s: %s" %
                                     (target, val))
                else:
                    installed_val = val
            elif len(installed_vals) == 1:
                installed_val = installed_vals[0]
            else:
                # check for partial matches
                installed_val = None
                for v in installed_vals:
                    if v.endswith(val + ".bed.gz") or v.endswith(val + ".bed"):
                        installed_val = v
                        break
                # handle date-stamped inputs
                if not installed_val:
                    installed_val = sorted(installed_vals, reverse=True)[0]
            data = tz.update_in(data, ["config", "algorithm", target], lambda x: installed_val)
    return data

# ## Sample and BAM read group naming


def _clean_metadata(data):
    batches = tz.get_in(("metadata", "batch"), data)
    # Ensure batches are strings and have no duplicates
    if batches:
        if isinstance(batches, (list, tuple)):
            batches = [_clean_characters(x) for x in sorted(list(set(batches)))]
        else:
            batches = _clean_characters(batches)
        data["metadata"]["batch"] = batches
    # If we have jointcalling, add a single batch if not present
    elif tz.get_in(["algorithm", "jointcaller"], data) or "gvcf" in tz.get_in(["algorithm", "tools_on"], data):
        if "metadata" not in data:
            data["metadata"] = {}
        data["metadata"]["batch"] = "%s-joint" % dd.get_sample_name(data)
    return data

def _clean_algorithm(data):
    """Clean algorithm keys, handling items that can be specified as lists or single items.
    """
    # convert single items to lists
    for key in ["variantcaller", "jointcaller", "svcaller"]:
        val = tz.get_in(["algorithm", key], data)
        if val:
            if not isinstance(val, (list, tuple)) and isinstance(val, six.string_types):
                val = [val]
            # check for cases like [false] or [None]
            if isinstance(val, (list, tuple)):
                if len(val) == 1 and not val[0] or (isinstance(val[0], six.string_types)
                                                    and val[0].lower() in ["none", "false"]):
                    val = False
            data["algorithm"][key] = val
    return data

def _organize_tools_on(data, is_cwl):
    """Ensure tools_on inputs match items specified elsewhere.
    """
    # want tools_on: [gvcf] if joint calling specified in CWL
    if is_cwl:
        if tz.get_in(["algorithm", "jointcaller"], data):
            val = tz.get_in(["algorithm", "tools_on"], data)
            if not val:
                val = []
            if not isinstance(val, (list, tuple)):
                val = [val]
            if "gvcf" not in val:
                val.append("gvcf")
            data["algorithm"]["tools_on"] = val
    return data

def _clean_background(data):
    """Clean up background specification, remaining back compatible.
    """
    allowed_keys = set(["variant", "cnv_reference"])
    val = tz.get_in(["algorithm", "background"], data)
    errors = []
    if val:
        out = {}
        # old style specification, single string for variant
        if isinstance(val, six.string_types):
            out["variant"] = _file_to_abs(val, [os.getcwd()])
        elif isinstance(val, dict):
            for k, v in val.items():
                if k in allowed_keys:
                    if isinstance(v, six.string_types):
                        out[k] = _file_to_abs(v, [os.getcwd()])
                    else:
                        assert isinstance(v, dict)
                        for ik, iv in v.items():
                            v[ik] = _file_to_abs(iv, [os.getcwd()])
                        out[k] = v
                else:
                    errors.append("Unexpected key: %s" % k)
        else:
            errors.append("Unexpected input: %s" % val)
        if errors:
            raise ValueError("Problematic algorithm background specification for %s:\n %s" %
                             (data["description"], "\n".join(errors)))
        out["cnv_reference"] = structural.standardize_cnv_reference({"config": data,
                                                                     "description": data["description"]})
        data["algorithm"]["background"] = out
    return data

def _clean_characters(x):
    """Clean problem characters in sample lane or descriptions.
    """
    if not isinstance(x, six.string_types):
        x = str(x)
    else:
        if not all(ord(char) < 128 for char in x):
            msg = "Found unicode character in input YAML (%s)" % (x)
            raise ValueError(repr(msg))
    for problem in [" ", ".", "/", "\\", "[", "]", "&", ";", "#", "+", ":", ")", "("]:
        x = x.replace(problem, "_")
    return x

def prep_rg_names(item, config, fc_name, fc_date):
    """Generate read group names from item inputs.
    """
    if fc_name and fc_date:
        lane_name = "%s_%s_%s" % (item["lane"], fc_date, fc_name)
    else:
        lane_name = item["description"]
    return {"rg": item["description"],
            "sample": item["description"],
            "lane": lane_name,
            "pl": (tz.get_in(["algorithm", "platform"], item)
                   or tz.get_in(["algorithm", "platform"], item, "illumina")).lower(),
            "lb": tz.get_in(["metadata", "library"], item),
            "pu": tz.get_in(["metadata", "platform_unit"], item) or lane_name}

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
        batches = tz.get_in(("metadata", "batch"), x)
        if batches:
            if not isinstance(batches, (list, tuple)):
                batches = [batches]
            for batch in batches:
                if batch in names:
                    dups.add(batch)
    if len(dups) > 0:
        raise ValueError("Batch names must be unique from sample descriptions.\n"
                         "Clashing batch names: %s" % sorted(list(dups)))

def _check_for_problem_somatic_batches(items, config):
    """Identify problem batch setups for somatic calling.

    We do not support multiple tumors in a single batch and VarDict(Java) does not
    handle pooled calling, only tumor/normal.
    """
    to_check = []
    for data in items:
        data = copy.deepcopy(data)
        data["config"] = config_utils.update_w_custom(config, data)
        to_check.append(data)
    data_by_batches = collections.defaultdict(list)
    for data in to_check:
        batches = dd.get_batches(data)
        if batches:
            for batch in batches:
                data_by_batches[batch].append(data)
    for batch, items in data_by_batches.items():
        if vcfutils.get_paired(items):
            vcfutils.check_paired_problems(items)
        elif len(items) > 1:
            vcs = vcfutils.get_somatic_variantcallers(items)
            if "vardict" in vcs:
                raise ValueError("VarDict does not support pooled non-tumor/normal calling, in batch %s: %s"
                                 % (batch, [dd.get_sample_name(data) for data in items]))
            elif "mutect" in vcs or "mutect2" in vcs:
                raise ValueError("MuTect and MuTect2 require a 'phenotype: tumor' sample for calling, "
                                 "in batch %s: %s"
                                 % (batch, [dd.get_sample_name(data) for data in items]))

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

def _check_for_degenerate_interesting_groups(items):
    """ Make sure interesting_groups specify existing metadata and that
    the interesting_group is not all of the same for all of the samples
    """
    igkey = ("algorithm", "bcbiornaseq", "interesting_groups")
    interesting_groups = tz.get_in(igkey, items[0], [])
    if isinstance(interesting_groups, str):
        interesting_groups = [interesting_groups]
    for group in interesting_groups:
        values = [tz.get_in(("metadata", group), x, None) for x in items]
        if all(x is None for x in values):
            raise ValueError("group %s is labelled as an interesting group, "
                             "but does not appear in the metadata." % group)
        if len(list(tz.unique(values))) == 1:
            raise ValueError("group %s is marked as an interesting group, "
                             "but all samples have the same value." % group)

TOPLEVEL_KEYS = set(["description", "analysis", "genome_build", "metadata", "algorithm",
                     "resources", "files", "vrn_file", "lane", "upload", "rgnames"])
ALGORITHM_KEYS = set(["bam_sort", "custom_trim", "kraken", "write_summary",
                      "merge_bamprep", "indelcaller", "effects", 
                      "svvalidate", "hlavalidate", "phasing", "validate",
                      "validate_regions", "validate_genome_build", "validate_method",
                      "clinical_reporting", "nomap_split_size", 
                      "nomap_split_targets", "background", "qc", "preseq",] +
                     # back compatibility
                     ["remove_lcr", "coverage_depth_max", "coverage_depth"] +
                     # from datadict.LOOKUPS
                     dd.get_algorithm_keys())
ALG_ALLOW_BOOLEANS = set(["merge_bamprep", "mark_duplicates", "remove_lcr",
                          "demultiplexed", "clinical_reporting", "transcriptome_align",
                          "fusion_mode", "assemble_transcripts", "trim_reads",
                          "quantify_genome_alignments", 
                          "recalibrate", "realign", "cwl_reporting", "save_diskspace"])
ALG_ALLOW_FALSE = set(["aligner", "align_split_size", "bam_clean", "bam_sort",
                       "effects", "phasing", "mixup_check", "indelcaller",
                       "variantcaller", "positional_umi", "maxcov_downsample", "preseq"])

ALG_DOC_URL = "https://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#algorithm-parameters"

def _check_algorithm_keys(item):
    """Check for unexpected keys in the algorithm section.

    Needs to be manually updated when introducing new keys, but avoids silent bugs
    with typos in key names.
    """
    problem_keys = [k for k in item["algorithm"].keys() if k not in ALGORITHM_KEYS]
    if len(problem_keys) > 0:
        raise ValueError("Unexpected configuration keyword in 'algorithm' section: %s\n"
                         "See configuration documentation for supported options:\n%s\n"
                         % (problem_keys, ALG_DOC_URL))

def _check_algorithm_values(item):
    """Check for misplaced inputs in the algorithms.

    - Identify incorrect boolean values where a choice is required.
    """
    problems = []
    for k, v in item.get("algorithm", {}).items():
        if v is True and k not in ALG_ALLOW_BOOLEANS:
            problems.append("%s set as true" % k)
        elif v is False and (k not in ALG_ALLOW_BOOLEANS and k not in ALG_ALLOW_FALSE):
            problems.append("%s set as false" % k)
    if len(problems) > 0:
        raise ValueError("Incorrect settings in 'algorithm' section for %s:\n%s"
                         "\nSee configuration documentation for supported options:\n%s\n"
                         % (item["description"], "\n".join(problems), ALG_DOC_URL))


def _check_toplevel_misplaced(item):
    """Check for algorithm keys accidentally placed at the top level.
    """
    problem_keys = [k for k in item.keys() if k in ALGORITHM_KEYS]
    if len(problem_keys) > 0:
        raise ValueError("Unexpected configuration keywords found in top level of %s: %s\n"
                         "This should be placed in the 'algorithm' section."
                         % (item["description"], problem_keys))
    problem_keys = [k for k in item.keys() if k not in TOPLEVEL_KEYS]
    if len(problem_keys) > 0:
        raise ValueError("Unexpected configuration keywords found in top level of %s: %s\n"
                         % (item["description"], problem_keys))

def _detect_fastq_format(in_file, MAX_RECORDS=1000):
    ranges = {"sanger": (33, 126),
              "solexa": (59, 126),
              "illumina_1.3+": (64, 126),
              "illumina_1.5+": (66, 126)}

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
            # if there is a short sequence, skip it
            if len(vals) < 20:
                continue
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
    fastq_extensions = ["fq.gz", "fastq.gz", ".fastq", ".fq"]

    for item in items:
        specified_format = item["algorithm"].get("quality_format", "standard").lower()
        if specified_format not in SAMPLE_FORMAT.values():
            raise ValueError("Quality format specified in the YAML file"
                             "is not supported. Supported values are %s."
                             % (SAMPLE_FORMAT.values()))

        fastq_file = next((f for f in item.get("files") or [] if f.endswith(tuple(fastq_extensions))), None)

        if fastq_file and specified_format and not objectstore.is_remote(fastq_file):
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
    allowed = set(list(alignment.TOOLS.keys()) + [None, False])
    if item["algorithm"].get("aligner") not in allowed:
        raise ValueError("Unexpected algorithm 'aligner' parameter: %s\n"
                         "Supported options: %s\n" %
                         (item["algorithm"].get("aligner"), sorted(list(allowed))))

def _check_variantcaller(item):
    """Ensure specified variantcaller is a valid choice.
    """
    allowed = set(list(genotype.get_variantcallers().keys()) + [None, False])
    vcs = item["algorithm"].get("variantcaller")
    if not isinstance(vcs, dict):
        vcs = {"variantcaller": vcs}
    for vc_set in vcs.values():
        if not isinstance(vc_set, (tuple, list)):
            vc_set = [vc_set]
        problem = [x for x in vc_set if x not in allowed]
        if len(problem) > 0:
            raise ValueError("Unexpected algorithm 'variantcaller' parameter: %s\n"
                             "Supported options: %s\n" % (problem, sorted(list(allowed))))
    # Ensure germline somatic calling only specified with tumor/normal samples
    if "germline" in vcs or "somatic" in vcs:
        paired = vcfutils.get_paired_phenotype(item)
        if not paired:
            raise ValueError("%s: somatic/germline calling in 'variantcaller' "
                             "but tumor/normal metadata phenotype not specified" % dd.get_sample_name(item))

def _check_svcaller(item):
    """Ensure the provide structural variant caller is valid.
    """
    allowed = set(reduce(operator.add, [list(d.keys()) for d in structural._CALLERS.values()]) + [None, False])
    svs = item["algorithm"].get("svcaller")
    if not isinstance(svs, (list, tuple)):
        svs = [svs]
    problem = [x for x in svs if x not in allowed]
    if len(problem) > 0:
        raise ValueError("Unexpected algorithm 'svcaller' parameters: %s\n"
                         "Supported options: %s\n" % (" ".join(["'%s'" % x for x in problem]),
                                                      sorted(list(allowed))))
    if "gatk-cnv" in svs and "cnvkit" in svs:
        raise ValueError("%s uses `gatk-cnv' and 'cnvkit', please use on one of these CNV callers" %
                         dd.get_sample_name(item))

def _get_as_list(item, k):
    out = item["algorithm"].get(k)
    if not out:
        out = []
    if not isinstance(out, (list, tuple)):
        out = [svs]
    return out

def _check_hetcaller(item):
    """Ensure upstream SV callers requires to heterogeneity analysis are available.
    """
    svs = _get_as_list(item, "svcaller")
    hets = _get_as_list(item, "hetcaller")
    if hets or any([x in svs for x in ["titancna", "purecn"]]):
        if not any([x in svs for x in ["cnvkit", "gatk-cnv"]]):
            raise ValueError("Heterogeneity caller used but need CNV calls. Add `gatk-cnv` "
                             "or `cnvkit` to `svcaller` in sample: %s" % item["description"])

def _check_jointcaller(data):
    """Ensure specified jointcaller is valid.
    """
    allowed = set(joint.get_callers() + [None, False])
    cs = data["algorithm"].get("jointcaller", [])
    if not isinstance(cs, (tuple, list)):
        cs = [cs]
    problem = [x for x in cs if x not in allowed]
    if len(problem) > 0:
        raise ValueError("Unexpected algorithm 'jointcaller' parameter: %s\n"
                         "Supported options: %s\n" % (problem, sorted(list(allowed), key=lambda x: x or "")))

def _check_indelcaller(data):
    c = data["algorithm"].get("indelcaller")
    if c and isinstance(c, (tuple, list)):
        raise ValueError("In sample %s, indelcaller specified as list. Can only be a single item: %s"
                         % (data["description"], str(c)))

def _check_hlacaller(data):
    supported_genomes = set(["hg38"])
    c = data["algorithm"].get("hlacaller")
    if c:
        if data["genome_build"] not in supported_genomes:
            raise ValueError("In sample %s, HLA caller specified but genome %s not in supported: %s" %
                             (data["description"], data["genome_build"], ", ".join(sorted(list(supported_genomes)))))

def _check_realign(data):
    """Check for realignment, which is not supported in GATK4
    """
    if "gatk4" not in data["algorithm"].get("tools_off", []) and not "gatk4" == data["algorithm"].get("tools_off"):
        if data["algorithm"].get("realign"):
            raise ValueError("In sample %s, realign specified but it is not supported for GATK4. "
                             "Realignment is generally not necessary for most variant callers." %
                             (dd.get_sample_name(data)))

def _check_trim(data):
    """Check for valid values for trim_reads.
    """
    trim = data["algorithm"].get("trim_reads")
    if trim:
        if trim == "fastp" and data["algorithm"].get("align_split_size") is not False:
            raise ValueError("In sample %s, `trim_reads: fastp` currently requires `align_split_size: false`" %
                             (dd.get_sample_name(data)))


def _check_sample_config(items, in_file, config):
    """Identify common problems in input sample configuration files.
    """
    logger.info("Checking sample YAML configuration: %s" % in_file)
    _check_quality_format(items)
    _check_for_duplicates(items, "lane")
    _check_for_duplicates(items, "description")
    _check_for_degenerate_interesting_groups(items)
    _check_for_batch_clashes(items)
    _check_for_problem_somatic_batches(items, config)
    _check_for_misplaced(items, "algorithm",
                         ["resources", "metadata", "analysis",
                          "description", "genome_build", "lane", "files"])

    [_check_toplevel_misplaced(x) for x in items]
    [_check_algorithm_keys(x) for x in items]
    [_check_algorithm_values(x) for x in items]
    [_check_aligner(x) for x in items]
    [_check_variantcaller(x) for x in items]
    [_check_svcaller(x) for x in items]
    [_check_hetcaller(x) for x in items]
    [_check_indelcaller(x) for x in items]
    [_check_jointcaller(x) for x in items]
    [_check_hlacaller(x) for x in items]
    [_check_realign(x) for x in items]
    [_check_trim(x) for x in items]

# ## Read bcbio_sample.yaml files

def _file_to_abs(x, dnames, makedir=False):
    """Make a file absolute using the supplied base directory choices.
    """
    if x is None or os.path.isabs(x):
        return x
    elif isinstance(x, six.string_types) and objectstore.is_remote(x):
        return x
    elif isinstance(x, six.string_types) and x.lower() == "none":
        return None
    else:
        for dname in dnames:
            if dname:
                normx = os.path.normpath(os.path.join(dname, x))
                if os.path.exists(normx):
                    return normx
                elif makedir:
                    utils.safe_makedir(normx)
                    return normx
        raise ValueError("Did not find input file %s in %s" % (x, dnames))

def _normalize_files(item, fc_dir=None):
    """Ensure the files argument is a list of absolute file names.
    Handles BAM, single and paired end fastq, as well as split inputs.
    """
    files = item.get("files")
    if files:
        if isinstance(files, six.string_types):
            files = [files]
        fastq_dir = flowcell.get_fastq_dir(fc_dir) if fc_dir else os.getcwd()
        files = [_file_to_abs(x, [os.getcwd(), fc_dir, fastq_dir]) for x in files]
        files = [x for x in files if x]
        _sanity_check_files(item, files)
        item["files"] = files
    return item

def _sanity_check_files(item, files):
    """Ensure input files correspond with supported approaches.

    Handles BAM, fastqs, plus split fastqs.
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
        if len(files) not in [1, 2] and item["analysis"].lower() != "scrna-seq":
            pair_types = set([len(xs) for xs in fastq.combine_pairs(files)])
            if len(pair_types) != 1 or pair_types.pop() not in [1, 2]:
                msg = "Expect either 1 (single end) or 2 (paired end) fastq inputs"
        if len(files) == 2 and files[0] == files[1]:
            msg = "Expect both fastq files to not be the same"
    if msg:
        raise ValueError("%s for %s: %s" % (msg, item.get("description", ""), files))

def validate_yaml(yaml_in, yaml_fn):
    """Check with yamllint the yaml syntaxes
    Looking for duplicate keys."""
    try:
        import yamllint.linter as linter
        from yamllint.config import YamlLintConfig
    except ImportError:
        return
    conf = """{"extends": "relaxed",
               "rules": {"trailing-spaces": {"level": "warning"},
                         "new-lines": {"level": "warning"},
                         "new-line-at-end-of-file": {"level": "warning"}}}"""
    if utils.file_exists(yaml_in):
        with open(yaml_in) as in_handle:
            yaml_in = in_handle.read()
    out = linter.run(yaml_in, YamlLintConfig(conf))

    for problem in out:
        msg = '%(fn)s:%(line)s:%(col)s: [%(level)s] %(msg)s' % {'fn': yaml_fn,
                                                                'line': problem.line,
                                                                'col': problem.column,
                                                                'level': problem.level,
                                                                'msg': problem.message}
        if problem.level == "error":
            raise ValueError(msg)

def _run_info_from_yaml(dirs, run_info_yaml, config, sample_names=None,
                        is_cwl=False, integrations=None):
    """Read run information from a passed YAML file.
    """
    validate_yaml(run_info_yaml, run_info_yaml)
    with open(run_info_yaml) as in_handle:
        loaded = yaml.safe_load(in_handle)
    fc_name, fc_date = None, None
    if dirs.get("flowcell"):
        try:
            fc_name, fc_date = flowcell.parse_dirname(dirs.get("flowcell"))
        except ValueError:
            pass
    global_config = {}
    global_vars = {}
    resources = {}
    integration_config = {}
    if isinstance(loaded, dict):
        global_config = copy.deepcopy(loaded)
        del global_config["details"]
        if "fc_name" in loaded:
            fc_name = loaded["fc_name"].replace(" ", "_")
        if "fc_date" in loaded:
            fc_date = str(loaded["fc_date"]).replace(" ", "_")
        global_vars = global_config.pop("globals", {})
        resources = global_config.pop("resources", {})
        for iname in ["arvados"]:
            integration_config[iname] = global_config.pop(iname, {})
        loaded = loaded["details"]
    if sample_names:
        loaded = [x for x in loaded if x["description"] in sample_names]

    if integrations:
        for iname, retriever in integrations.items():
            if iname in config:
                config[iname] = retriever.set_cache(config[iname])
                loaded = retriever.add_remotes(loaded, config[iname])

    run_details = []
    for i, item in enumerate(loaded):
        item = _normalize_files(item, dirs.get("flowcell"))
        if "lane" not in item:
            item["lane"] = str(i + 1)
        item["lane"] = _clean_characters(item["lane"])
        if "description" not in item:
            if _item_is_bam(item):
                item["description"] = get_sample_name(item["files"][0])
            else:
                raise ValueError("No `description` sample name provided for input #%s" % (i + 1))
        description = _clean_characters(item["description"])
        item["description"] = description
        # make names R safe if we are likely to use R downstream
        if item["analysis"].lower() in R_DOWNSTREAM_ANALYSIS:
            if description[0].isdigit():
                valid = "X" + description
                logger.info("%s is not a valid R name, converting to %s." % (description, valid))
                item["description"] = valid
        if "upload" not in item and not is_cwl:
            upload = global_config.get("upload", {})
            # Handle specifying a local directory directly in upload
            if isinstance(upload, six.string_types):
                upload = {"dir": upload}
            if not upload:
                upload["dir"] = "../final"
            if fc_name:
                upload["fc_name"] = fc_name
            if fc_date:
                upload["fc_date"] = fc_date
            upload["run_id"] = ""
            if upload.get("dir"):
                upload["dir"] = _file_to_abs(upload["dir"], [dirs.get("work")], makedir=True)
            item["upload"] = upload
        item["algorithm"] = _replace_global_vars(item["algorithm"], global_vars)
        item["algorithm"] = genome.abs_file_paths(item["algorithm"],
                                                  ignore_keys=ALGORITHM_NOPATH_KEYS,
                                                  fileonly_keys=ALGORITHM_FILEONLY_KEYS,
                                                  do_download=all(not x for x in integrations.values()))
        item["genome_build"] = str(item.get("genome_build", ""))
        item["algorithm"] = _add_algorithm_defaults(item["algorithm"], item.get("analysis", ""), is_cwl)
        item["metadata"] = add_metadata_defaults(item.get("metadata", {}))
        item["rgnames"] = prep_rg_names(item, config, fc_name, fc_date)
        if item.get("files"):
            item["files"] = [genome.abs_file_paths(f, do_download=all(not x for x in integrations.values()))
                             for f in item["files"]]
        elif "files" in item:
            del item["files"]
        if item.get("vrn_file") and isinstance(item["vrn_file"], six.string_types):
            item["vrn_file"] = genome.abs_file_paths(item["vrn_file"],
                                                     do_download=all(not x for x in integrations.values()))
            if os.path.isfile(item["vrn_file"]):
                # Try to prepare in place (or use ready to go inputs)
                try:
                    item["vrn_file"] = vcfutils.bgzip_and_index(item["vrn_file"], config,
                                                                remove_orig=False)
                # In case of permission errors, fix in inputs directory
                except IOError:
                    inputs_dir = utils.safe_makedir(os.path.join(dirs.get("work", os.getcwd()), "inputs",
                                                                 item["description"]))
                    item["vrn_file"] = vcfutils.bgzip_and_index(item["vrn_file"], config,
                                                                remove_orig=False, out_dir=inputs_dir)
            if not tz.get_in(("metadata", "batch"), item) and tz.get_in(["algorithm", "validate"], item):
                raise ValueError("%s: Please specify a metadata batch for variant file (vrn_file) input.\n" %
                                 (item["description"]) +
                                 "Batching with a standard sample provides callable regions for validation.")
        item = _clean_metadata(item)
        item = _clean_algorithm(item)
        item = _organize_tools_on(item, is_cwl)
        item = _clean_background(item)
        # Add any global resource specifications
        if "resources" not in item:
            item["resources"] = {}
        for prog, pkvs in resources.items():
            if prog not in item["resources"]:
                item["resources"][prog] = {}
            if pkvs is not None:
                for key, val in pkvs.items():
                    item["resources"][prog][key] = val
        for iname, ivals in integration_config.items():
            if ivals:
                if iname not in item:
                    item[iname] = {}
                for k, v in ivals.items():
                    item[iname][k] = v

        run_details.append(item)
    _check_sample_config(run_details, run_info_yaml, config)
    return run_details

def _item_is_bam(item):
    files = item.get("files", [])
    return len(files) == 1 and files[0].endswith(".bam")

def add_metadata_defaults(md):
    """Central location for defaults for algorithm inputs.
    """
    defaults = {"batch": None,
                "phenotype": ""}
    for k, v in defaults.items():
        if k not in md:
            md[k] = v
    return md

def _get_nomap_split_targets(analysis, is_cwl):
    """Chromosome splitting logic based on run type.

    RNA-seq -- aim for smaller chunks (half chromosomes) to avoid memory issues
    CWL -- aim for larger chunks to allow batching and multicore
    old style -- larger number of chunks for better parallelization
    """
    if analysis.lower().find("rna-seq") >= 0:
        return 50
    elif is_cwl:
        return 20
    else:
        return 200

def _add_algorithm_defaults(algorithm, analysis, is_cwl):
    """Central location specifying defaults for algorithm inputs.

    Converts allowed multiple inputs into lists if specified as a single item.
    Converts required single items into string if specified as a list
    """
    if not algorithm:
        algorithm = {}
    defaults = {"archive": None,
                "tools_off": [],
                "tools_on": [],
                "qc": [],
                "trim_reads": False,
                "adapters": [],
                "effects": "snpeff",
                "quality_format": "standard",
                "expression_caller": ["salmon"] if analysis.lower().find("rna-seq") >= 0 else None,
                "align_split_size": None,
                "bam_clean": False,
                "nomap_split_size": 250,
                "nomap_split_targets": _get_nomap_split_targets(analysis, is_cwl),
                "mark_duplicates": False if not algorithm.get("aligner") else True,
                "coverage_interval": None,
                "min_allele_fraction": 10.0,
                "recalibrate": False,
                "realign": False,
                "ensemble": None,
                "exclude_regions": [],
                "variant_regions": None,
                "svcaller": [],
                "svvalidate": None,
                "svprioritize": None,
                "validate": None,
                "validate_regions": None,
                "vcfanno": []}
    convert_to_list = set(["tools_off", "tools_on", "hetcaller", "variantcaller", "svcaller", "qc", "disambiguate",
                           "vcfanno", "adapters", "custom_trim", "exclude_regions"])
    convert_to_single = set(["hlacaller", "indelcaller", "validate_method"])
    for k, v in defaults.items():
        if k not in algorithm:
            algorithm[k] = v
    for k, v in algorithm.items():
        if k in convert_to_list:
            if v and not isinstance(v, (list, tuple)) and not isinstance(v, dict):
                algorithm[k] = [v]
            # ensure dictionary specified inputs get converted into individual lists
            elif v and not isinstance(v, (list, tuple)) and isinstance(v, dict):
                new = {}
                for innerk, innerv in v.items():
                    if innerv and not isinstance(innerv, (list, tuple)) and not isinstance(innerv, dict):
                        innerv = [innerv]
                    new[innerk] = innerv
                algorithm[k] = new
            elif v is None:
                algorithm[k] = []
        elif k in convert_to_single:
            if v and not isinstance(v, six.string_types):
                if isinstance(v, (list, tuple)) and len(v) == 1:
                    algorithm[k] = v[0]
                else:
                    raise ValueError("Unexpected input in sample YAML; need a single item for %s: %s" % (k, v))
    return algorithm

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
        for k, v in xs.items():
            if isinstance(v, six.string_types) and v in global_vars:
                v = global_vars[v]
            final[k] = v
        return final
    else:
        return xs

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

def prep_system(run_info_yaml, bcbio_system=None):
    """Prepare system configuration information from an input configuration file.

    This does the work of parsing the system input file and setting up directories
    for use in 'organize'.
    """
    work_dir = os.getcwd()
    config, config_file = config_utils.load_system_config(bcbio_system, work_dir)
    dirs = setup_directories(work_dir, os.path.normpath(os.path.dirname(os.path.dirname(run_info_yaml))),
                             config, config_file)
    return [dirs, config, run_info_yaml]
