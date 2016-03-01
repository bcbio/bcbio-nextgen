"""Next-gen variant detection and evaluation with GATK and SnpEff.
"""
import os
import toolz as tz

from bcbio import utils
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.variation.genotype import variant_filtration, get_variantcaller
from bcbio.variation import effects, genotype, germline, prioritize
from bcbio.variation import multi as vmulti

# ## Genotyping

def postprocess_variants(items):
    """Provide post-processing of variant calls: filtering and effects annotation.
    """
    data = _get_batch_representative(items, "vrn_file")
    cur_name = "%s, %s" % (dd.get_sample_name(data), get_variantcaller(data))
    logger.info("Finalizing variant calls: %s" % cur_name)
    orig_vrn_file = data.get("vrn_file")
    data = _symlink_to_workdir(data, ["vrn_file"])
    data = _symlink_to_workdir(data, ["config", "algorithm", "variant_regions"])
    if data.get("align_bam") and data.get("vrn_file"):
        logger.info("Calculating variation effects for %s" % cur_name)
        ann_vrn_file, vrn_stats = effects.add_to_vcf(data["vrn_file"], data)
        if ann_vrn_file:
            data["vrn_file"] = ann_vrn_file
        if vrn_stats:
            data["vrn_stats"] = vrn_stats
        logger.info("Filtering for %s" % cur_name)
        orig_items = _get_orig_items(items)
        data["vrn_file"] = variant_filtration(data["vrn_file"], dd.get_ref_file(data),
                                              tz.get_in(("genome_resources", "variation"), data, {}),
                                              data, orig_items)
        logger.info("Prioritization for %s" % cur_name)
        data["vrn_file"] = prioritize.handle_vcf_calls(data["vrn_file"], data, orig_items)
        logger.info("Germline extraction for %s" % cur_name)
        data = germline.extract(data, orig_items)
    if orig_vrn_file and os.path.samefile(data["vrn_file"], orig_vrn_file):
        data["vrn_file"] = orig_vrn_file
    return [[data]]

def _get_orig_items(data):
    """Retrieve original items in a batch, handling CWL and standard cases.
    """
    if isinstance(data, dict):
        if tz.get_in(["metadata", "batch"], data):
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
        variantcaller = genotype.get_variantcaller(data)
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
        return items
    else:
        vals = set([])
        out = []
        for data in items:
            if key in data:
                vals.add(data[key])
                out.append(data)
        if len(vals) != 1:
            raise ValueError("Incorrect values for %s: %s" % (key, list(vals)))
        return out[0]
