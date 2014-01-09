"""Perform validation of final calls against known reference materials.

Automates the process of checking pipeline results against known valid calls
to identify discordant variants. This provides a baseline for ensuring the
validity of pipeline updates and algorithm changes.
"""
import csv
import os

import yaml

from bcbio import utils
from bcbio.bam import callable
from bcbio.pipeline import config_utils
from bcbio.pipeline import genome
from bcbio.provenance import do

# ## Individual sample comparisons

def _has_validate(data):
    return data.get("vrn_file") and "validate" in data["config"]["algorithm"]

def normalize_input_path(x, data):
    """Normalize path for input files, handling relative paths.
    Looks for non-absolute paths in local and fastq directories
    """
    if x is None:
        return None
    elif os.path.isabs(x):
        return os.path.normpath(x)
    else:
        for d in [data["dirs"].get("fastq"), data["dirs"].get("work")]:
            if d:
                cur_x = os.path.normpath(os.path.join(d, x))
                if os.path.exists(cur_x):
                    return cur_x
        raise IOError("Could not find validation file %s" % x)

def compare_to_rm(data):
    """Compare final variant calls against reference materials of known calls.
    """
    if _has_validate(data):
        if isinstance(data["vrn_file"], (list, tuple)):
            vrn_file = [os.path.abspath(x) for x in data["vrn_file"]]
        else:
            vrn_file = os.path.abspath(data["vrn_file"])
        rm_file = normalize_input_path(data["config"]["algorithm"]["validate"], data)
        rm_interval_file = normalize_input_path(data["config"]["algorithm"].get("validate_regions"), data)
        rm_genome = data["config"]["algorithm"].get("validate_genome_build")
        sample = data["name"][-1].replace(" ", "_")
        caller = data["config"]["algorithm"].get("variantcaller")
        if not caller:
            caller = "precalled"
        base_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "validate", sample, caller))
        val_config_file = _create_validate_config_file(vrn_file, rm_file, rm_interval_file,
                                                       rm_genome, base_dir, data)
        work_dir = os.path.join(base_dir, "work")
        out = {"summary": os.path.join(work_dir, "validate-summary.csv"),
               "grading": os.path.join(work_dir, "validate-grading.yaml"),
               "concordant": os.path.join(work_dir, "%s-ref-eval-concordance.vcf" % sample),
               "discordant": os.path.join(work_dir, "%s-eval-ref-discordance-annotate.vcf" % sample)}
        if not utils.file_exists(out["concordant"]) or not utils.file_exists(out["grading"]):
            bcbio_variation_comparison(val_config_file, base_dir, data)
        data["validate"] = out
    return [[data]]

def bcbio_variation_comparison(config_file, base_dir, data):
    """Run a variant comparison using the bcbio.variation toolkit, given an input configuration.
    """
    tmp_dir = utils.safe_makedir(os.path.join(base_dir, "tmp"))
    bv_jar = config_utils.get_jar("bcbio.variation",
                                  config_utils.get_program("bcbio_variation",
                                                           data["config"], "dir"))
    resources = config_utils.get_resources("bcbio_variation", data["config"])
    jvm_opts = resources.get("jvm_opts", ["-Xms750m", "-Xmx2g"])
    java_args = ["-Djava.io.tmpdir=%s" % tmp_dir]
    cmd = ["java"] + jvm_opts + java_args + ["-jar", bv_jar, "variant-compare", config_file]
    do.run(cmd, "Comparing variant calls using bcbio.variation", data)

def _create_validate_config_file(vrn_file, rm_file, rm_interval_file, rm_genome,
                                 base_dir, data):
    config_dir = utils.safe_makedir(os.path.join(base_dir, "config"))
    config_file = os.path.join(config_dir, "validate.yaml")
    with open(config_file, "w") as out_handle:
        out = _create_validate_config(vrn_file, rm_file, rm_interval_file, rm_genome,
                                      base_dir, data)
        yaml.dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return config_file

def _create_validate_config(vrn_file, rm_file, rm_interval_file, rm_genome,
                            base_dir, data):
    """Create a bcbio.variation configuration input for validation.
    """
    if rm_genome:
        rm_genome = genome.get_refs(rm_genome, None, data["dirs"]["galaxy"])[-1]
        if rm_genome != data["sam_ref"]:
            eval_genome = data["sam_ref"]
        else:
            eval_genome = None
    else:
        eval_genome = None
        rm_genome = data["sam_ref"]
    ref_call = {"file": str(rm_file), "name": "ref", "type": "grading-ref",
                "preclean": True, "prep": True, "remove-refcalls": True}
    a_intervals = get_analysis_intervals(data)
    if rm_interval_file:
        ref_call["intervals"] = rm_interval_file
    eval_call = {"file": vrn_file, "name": "eval", "remove-refcalls": True}
    if eval_genome:
        eval_call["ref"] = eval_genome
        eval_call["preclean"] = True
        eval_call["prep"] = True
    if a_intervals and eval_genome:
        eval_call["intervals"] = os.path.abspath(a_intervals)
    exp = {"sample": data["name"][-1],
           "ref": rm_genome,
           "approach": "grade",
           "calls": [ref_call, eval_call]}
    if a_intervals and not eval_genome:
        exp["intervals"] = os.path.abspath(a_intervals)
    if data.get("callable_bam") and not eval_genome:
        exp["align"] = data["callable_bam"]
    return {"dir": {"base": base_dir, "out": "work", "prep": "work/prep"},
            "experiments": [exp]}

def get_analysis_intervals(data):
    """Retrieve analysis regions for the current variant calling pipeline.
    """
    if data.get("ensemble_bed"):
        return data["ensemble_bed"]
    elif data.get("callable_bam"):
        return callable.sample_callable_bed(data["callable_bam"], data["sam_ref"], data["config"])
    else:
        for key in ["callable_regions", "variant_regions"]:
            intervals = data["config"]["algorithm"].get(key)
            if intervals:
                return intervals

# ## Summarize comparisons

def _flatten_grading(stats):
    vtypes = ["snp", "indel"]
    cat = "concordant"
    for vtype in vtypes:
        yield vtype, cat, stats[cat][cat].get(vtype, 0)
    for vtype in vtypes:
        for vclass, vitems in sorted(stats["discordant"].get(vtype, {}).iteritems()):
            for vreason, val in sorted(vitems.iteritems()):
                yield vtype, "discordant-%s-%s" % (vclass, vreason), val
            yield vtype, "discordant-%s-total" % vclass, sum(vitems.itervalues())

def _has_grading_info(samples):
    for data in (x[0] for x in samples):
        for variant in data.get("variants", []):
            if "validate" in variant:
                return True
    return False

def summarize_grading(samples):
    """Provide summaries of grading results across all samples.
    """
    if not _has_grading_info(samples):
        return samples
    out = []
    out_csv = os.path.join(samples[0][0]["dirs"]["work"],
                           "grading-summary.csv")
    with open(out_csv, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["sample", "caller", "variant.type", "category", "value"])
        for data in (x[0] for x in samples):
            for variant in data.get("variants", []):
                if variant.get("validate"):
                    variant["validate"]["grading_summary"] = out_csv
                    with open(variant["validate"]["grading"]) as in_handle:
                        grade_stats = yaml.load(in_handle)
                    for sample_stats in grade_stats:
                        sample = sample_stats["sample"]
                        for vtype, cat, val in _flatten_grading(sample_stats):
                            writer.writerow([sample, variant.get("variantcaller", ""),
                                             vtype, cat, val])
            out.append([data])
    return out
