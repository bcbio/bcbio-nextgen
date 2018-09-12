"""Next-gen variant detection and evaluation with GATK and SnpEff.
"""
import collections
import os
import toolz as tz

from bcbio import utils
from bcbio.cwl import cwlutils
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.variation.genotype import variant_filtration, get_variantcaller
from bcbio.variation import (annotation, damage, effects, genotype, germline, population, prioritize,
                             validate, vcfutils)
from bcbio.variation import multi as vmulti

# ## CWL summarization

def summarize_vc(items):
    """CWL target: summarize variant calls and validation for multiple samples.
    """
    items = [utils.to_single_data(x) for x in utils.flatten(items)]
    items = [_normalize_vc_input(x) for x in items]
    items = validate.summarize_grading(items)
    items = [utils.to_single_data(x) for x in items]
    out = {"validate": validate.combine_validations(items),
           "variants": {"calls": [], "gvcf": [], "samples": []}}
    added = set([])
    variants_by_sample = collections.defaultdict(list)
    sample_order = []
    for data in items:
        batch_samples = data.get("batch_samples", [dd.get_sample_name(data)])
        for s in batch_samples:
            if s not in sample_order:
                sample_order.append(s)
        if data.get("vrn_file"):
            # Only get batches if we're actually doing variantcalling in bcbio
            # otherwise we'll be using the original files
            names = dd.get_batches(data) if dd.get_variantcaller(data) else None
            if not names:
                names = [dd.get_sample_name(data)]
            batch_name = names[0]
            if data.get("vrn_file_joint") is not None:
                to_add = [("vrn_file", "gvcf", dd.get_sample_name(data)),
                          ("vrn_file_joint", "calls", batch_name)]
            else:
                to_add = [("vrn_file", "calls", batch_name)]
            for vrn_key, out_key, name in to_add:
                cur_name = "%s-%s" % (name, dd.get_variantcaller(data))
                out_file = os.path.join(utils.safe_makedir(os.path.join(dd.get_work_dir(data),
                                                                        "variants", out_key)),
                                        "%s.vcf.gz" % cur_name)
                for s in batch_samples:
                    variants_by_sample[s].append(out_file)
                if cur_name not in added:
                    added.add(cur_name)
                    # Ideally could symlink here but doesn't appear to work with
                    # Docker container runs on Toil where PATHs don't get remapped
                    utils.copy_plus(os.path.realpath(data[vrn_key]), out_file)
                    vcfutils.bgzip_and_index(out_file, data["config"])
                    out["variants"][out_key].append(out_file)
    for sample in sample_order:
        out["variants"]["samples"].append(variants_by_sample[sample])
    return [out]

def _normalize_vc_input(data):
    """Normalize different types of variant calling inputs.

    Handles standard and ensemble inputs.
    """
    if data.get("ensemble"):
        for k in ["batch_samples", "validate", "vrn_file"]:
            data[k] = data["ensemble"][k]
        data["config"]["algorithm"]["variantcaller"] = "ensemble"
        data["metadata"] = {"batch": data["ensemble"]["batch_id"]}
    return data

# ## Genotyping

def postprocess_variants(items):
    """Provide post-processing of variant calls: filtering and effects annotation.
    """
    vrn_key = "vrn_file"
    if not isinstance(items, dict):
        items = [utils.to_single_data(x) for x in items]
        if "vrn_file_joint" in items[0]:
            vrn_key = "vrn_file_joint"
    data, items = _get_batch_representative(items, vrn_key)
    items = cwlutils.unpack_tarballs(items, data)
    data = cwlutils.unpack_tarballs(data, data)
    cur_name = "%s, %s" % (dd.get_sample_name(data), get_variantcaller(data, require_bam=False))
    logger.info("Finalizing variant calls: %s" % cur_name)
    orig_vrn_file = data.get(vrn_key)
    data = _symlink_to_workdir(data, [vrn_key])
    data = _symlink_to_workdir(data, ["config", "algorithm", "variant_regions"])
    if data.get(vrn_key):
        logger.info("Calculating variation effects for %s" % cur_name)
        ann_vrn_file, vrn_stats = effects.add_to_vcf(data[vrn_key], data)
        if ann_vrn_file:
            data[vrn_key] = ann_vrn_file
        if vrn_stats:
            data["vrn_stats"] = vrn_stats
        orig_items = _get_orig_items(items)
        logger.info("Annotate VCF file: %s" % cur_name)
        data[vrn_key] = annotation.finalize_vcf(data[vrn_key], get_variantcaller(data, require_bam=False),
                                                orig_items)
        if cwlutils.is_cwl_run(data):
            logger.info("Annotate with population level variation data")
            ann_file = population.run_vcfanno(data[vrn_key], data)
            if ann_file:
                data[vrn_key] = ann_file
        logger.info("Filtering for %s" % cur_name)
        data[vrn_key] = variant_filtration(data[vrn_key], dd.get_ref_file(data),
                                           tz.get_in(("genome_resources", "variation"), data, {}),
                                           data, orig_items)
        logger.info("Prioritization for %s" % cur_name)
        prio_vrn_file = prioritize.handle_vcf_calls(data[vrn_key], data, orig_items)
        if prio_vrn_file != data[vrn_key]:
            data[vrn_key] = prio_vrn_file
            logger.info("Germline extraction for %s" % cur_name)
            data = germline.extract(data, orig_items)

        if dd.get_align_bam(data):
            data = damage.run_filter(data[vrn_key], dd.get_align_bam(data), dd.get_ref_file(data),
                                     data, orig_items)
    if orig_vrn_file and os.path.samefile(data[vrn_key], orig_vrn_file):
        data[vrn_key] = orig_vrn_file
    return [[data]]

def _get_orig_items(data):
    """Retrieve original items in a batch, handling CWL and standard cases.
    """
    if isinstance(data, dict):
        if dd.get_align_bam(data) and tz.get_in(["metadata", "batch"], data) and "group_orig" in data:
            return vmulti.get_orig_items(data)
        else:
            return [data]
    else:
        return data

def _symlink_to_workdir(data, key):
    """For CWL support, symlink files into a working directory if in read-only imports.
    """
    orig_file = tz.get_in(key, data)
    if orig_file and not orig_file.startswith(dd.get_work_dir(data)):
        variantcaller = genotype.get_variantcaller(data, require_bam=False)
        if not variantcaller:
            variantcaller = "precalled"
        out_file = os.path.join(dd.get_work_dir(data), variantcaller, os.path.basename(orig_file))
        utils.safe_makedir(os.path.dirname(out_file))
        utils.symlink_plus(orig_file, out_file)
        data = tz.update_in(data, key, lambda x: out_file)
    return data

def _get_batch_representative(items, key):
    """Retrieve a representative data item from a batch.

    Handles standard bcbio cases (a single data item) and CWL cases with
    batches that have a consistent variant file.
    """
    if isinstance(items, dict):
        return items, items
    else:
        vals = set([])
        out = []
        for data in items:
            if key in data:
                vals.add(data[key])
                out.append(data)
        if len(vals) != 1:
            raise ValueError("Incorrect values for %s: %s" % (key, list(vals)))
        return out[0], items
