"""Organize samples for coordinated multi-sample processing.

Handles grouping of related families or batches to go through variant
calling simultaneously.
"""
import collections
import copy
import os

from bcbio import utils
from bcbio.variation import vcfutils

def group_by_batch(items):
    """Group a set of sample items by batch (or singleton) name.

    Items in multiple batches cause two batches to be merged together.
    """
    out = collections.defaultdict(list)
    batch_groups = _get_representative_batch(_merge_batches(_find_all_groups(items)))
    for data in items:
        batch = utils.get_in(data, ("metadata", "batch"), data["description"])
        if isinstance(batch, (list, tuple)):
            batch = batch[0]
        batch = batch_groups[batch]
        out[batch].append(data)
    return dict(out)

def _find_all_groups(items):
    """Find all groups
    """
    all_groups = []
    for data in items:
        batches = utils.get_in(data, ("metadata", "batch"), data["description"])
        if not isinstance(batches, (list, tuple)):
            batches = [batches]
        all_groups.append(batches)
    return all_groups

def _merge_batches(all_groups):
    """Merge batches with overlapping groups. Uses merge approach from:

    http://stackoverflow.com/a/4842897/252589
    """
    merged = []
    while len(all_groups) > 0:
        first, rest = all_groups[0], all_groups[1:]
        first = set(first)
        lf = -1
        while len(first) > lf:
            lf = len(first)

            rest2 = []
            for r in rest:
                if len(first.intersection(set(r))) > 0:
                    first |= set(r)
                else:
                    rest2.append(r)
            rest = rest2
        merged.append(first)
        all_groups = rest
    return merged

def _get_representative_batch(merged):
    """Prepare dictionary matching batch items to a representative within a group.
    """
    out = {}
    for mgroup in merged:
        mgroup = sorted(list(mgroup))
        for x in mgroup:
            out[x] = mgroup[0]
    return out

def group_batches(xs):
    """Group samples into batches for simultaneous variant calling.

    Identify all samples to call together: those in the same batch
    and variant caller.
    Pull together all BAM files from this batch and process together,
    Provide details to pull these finalized files back into individual
    expected files.
    """
    singles = []
    batch_groups = collections.defaultdict(list)
    for args in xs:
        assert len(args) == 1
        data = args[0]
        batch = utils.get_in(data, ("metadata", "batch"))
        caller = data["config"]["algorithm"].get("variantcaller", "gatk")
        region = tuple(data["region"]) if "region" in data else ()
        if batch is not None:
            batches = batch if isinstance(batch, (list, tuple)) else [batch]
            for b in batches:
                batch_groups[(b, region, caller)].append(data)
        else:
            singles.append(data)
    batches = []
    for batch, items in batch_groups.iteritems():
        batch_data = copy.deepcopy(_pick_lead_item(items))
        batch_data["work_bam"] = [x["work_bam"] for x in items]
        batch_data["group_orig"] = items
        batch_data["group"] = batch
        batches.append(batch_data)
    return singles + batches

def _pick_lead_item(items):
    """Pick single representative sample for batch calling to attach calls to.

    For cancer samples, attach to tumor.
    """
    if vcfutils.is_paired_analysis([x["work_bam"] for x in items], items):
        for data in items:
            if vcfutils.get_paired_phenotype(data) == "tumor":
                return data
        raise ValueError("Did not find tumor sample in paired tumor/normal calling")
    else:
        return items[0]

def split_variants_by_sample(data):
    """Split a multi-sample call file into inputs for individual samples.

    For tumor/normal paired analyses, do not split the final file and attach
    it to the tumor input.
    """
    # not split, do nothing
    if "group_orig" not in data:
        return [[data]]
    # cancer tumor/normal
    elif vcfutils.get_paired_phenotype(data):
        out = []
        for i, sub_data in enumerate(data["group_orig"]):
            if vcfutils.get_paired_phenotype(sub_data) == "tumor":
                sub_data["vrn_file"] = data["vrn_file"]
            out.append([sub_data])
        return out
    # population or single sample
    else:
        out = []
        for sub_data in data["group_orig"]:
            sub_vrn_file = data["vrn_file"].replace(data["group"][0] + "-", sub_data["name"][-1] + "-")
            if len(vcfutils.get_samples(data["vrn_file"])) > 1:
                vcfutils.select_sample(data["vrn_file"], sub_data["name"][-1], sub_vrn_file, data["config"])
            elif not os.path.exists(sub_vrn_file):
                utils.symlink_plus(data["vrn_file"], sub_vrn_file)
            sub_data["vrn_file"] = sub_vrn_file
            out.append([sub_data])
        return out
