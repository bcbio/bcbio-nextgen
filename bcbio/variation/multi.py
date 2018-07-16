"""Organize samples for coordinated multi-sample processing.

Handles grouping of related families or batches to go through variant
calling simultaneously.
"""
import collections

import toolz as tz

from bcbio import utils
from bcbio.pipeline import datadict as dd
from bcbio.variation import vcfutils

# ## Group batches to process together

def group_by_batch(items, require_bam=True):
    """Group a set of sample items by batch (or singleton) name.

    Items in multiple batches cause two batches to be merged together.
    """
    out = collections.defaultdict(list)
    batch_groups = _get_representative_batch(_merge_batches(_find_all_groups(items, require_bam)))
    for data in items:
        batches = _get_batches(data, require_bam)
        # take first batch as representative
        batch = batch_groups[batches[0]]
        out[batch].append(utils.deepish_copy(data))
    return dict(out)

def bam_needs_processing(data):
    """Check if a work input needs processing for parallelization.
    """
    return ((data.get("work_bam") or data.get("align_bam")) and
            (any(tz.get_in(["config", "algorithm", x], data) for x in
                 ["variantcaller", "mark_duplicates", "recalibrate", "realign", "svcaller",
                  "jointcaller", "variant_regions"])
             or any(k in data for k in ["cwl_keys", "output_cwl_keys"])))

def get_batch_for_key(data):
    """Retrieve batch information useful as a unique key for the sample.
    """
    batches = _get_batches(data, require_bam=False)
    if len(batches) == 1:
        return batches[0]
    else:
        return tuple(batches)

def _get_batches(data, require_bam=True):
    if bam_needs_processing(data) or not require_bam:
        batches = dd.get_batch(data) or dd.get_sample_name(data)
    else:
        batches = dd.get_sample_name(data)
    if not isinstance(batches, (list, tuple)):
        batches = [batches]
    return batches

def _find_all_groups(items, require_bam=True):
    """Find all groups
    """
    all_groups = []
    for data in items:
        batches = _get_batches(data, require_bam)
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

def _list_to_tuple(xs):
    if isinstance(xs, (list, tuple)):
        return tuple([_list_to_tuple(x) for x in xs])
    else:
        return xs

def _group_batches_shared(xs, caller_batch_fn, prep_data_fn):
    """Shared functionality for grouping by batches for variant calling and joint calling.
    """
    singles = []
    batch_groups = collections.defaultdict(list)
    for args in xs:
        data = utils.to_single_data(args)
        caller, batch = caller_batch_fn(data)
        region = _list_to_tuple(data["region"]) if "region" in data else ()
        if batch is not None:
            batches = batch if isinstance(batch, (list, tuple)) else [batch]
            for b in batches:
                batch_groups[(b, region, caller)].append(utils.deepish_copy(data))
        else:
            data = prep_data_fn(data, [data])
            singles.append(data)
    batches = []
    for batch, items in batch_groups.items():
        batch_data = utils.deepish_copy(_pick_lead_item(items))
        # For nested primary batches, split permanently by batch
        if tz.get_in(["metadata", "batch"], batch_data):
            batch_name = batch[0]
            batch_data["metadata"]["batch"] = batch_name
        batch_data = prep_data_fn(batch_data, items)
        batch_data["group_orig"] = _collapse_subitems(batch_data, items)
        batch_data["group"] = batch
        batches.append(batch_data)
    return singles + batches

def group_batches(xs):
    """Group samples into batches for simultaneous variant calling.

    Identify all samples to call together: those in the same batch
    and variant caller.
    Pull together all BAM files from this batch and process together,
    Provide details to pull these finalized files back into individual
    expected files.
    Only batches files if joint calling not specified.
    """
    def _caller_batches(data):
        caller = tz.get_in(("config", "algorithm", "variantcaller"), data)
        jointcaller = tz.get_in(("config", "algorithm", "jointcaller"), data)
        batch = tz.get_in(("metadata", "batch"), data) if not jointcaller else None
        return caller, batch
    def _prep_data(data, items):
        data["region_bams"] = [x["region_bams"] for x in items]
        return data
    return _group_batches_shared(xs, _caller_batches, _prep_data)

def group_batches_joint(samples):
    """Perform grouping by batches for joint calling/squaring off.
    """
    def _caller_batches(data):
        jointcaller = tz.get_in(("config", "algorithm", "jointcaller"), data)
        batch = tz.get_in(("metadata", "batch"), data) if jointcaller else None
        return jointcaller, batch
    def _prep_data(data, items):
        for r in ["callable_regions", "variant_regions"]:
            data[r] = list(set(filter(lambda x: x is not None,
                                      [tz.get_in(("config", "algorithm", r), d) for d in items])))
        data["work_bams"] = [dd.get_align_bam(x) or dd.get_work_bam(x) for x in items]
        data["vrn_files"] = [x["vrn_file"] for x in items]
        return data
    return _group_batches_shared(samples, _caller_batches, _prep_data)

# ## Collapse and uncollapse groups to save memory

def _collapse_subitems(base, items):
    """Collapse full data representations relative to a standard base.
    """
    out = []
    for d in items:
        newd = _diff_dict(base, d)
        out.append(newd)
    return out

def _diff_dict(orig, new):
    """Diff a nested dictionary, returning only key/values that differ.
    """
    final = {}
    for k, v in new.items():
        if isinstance(v, dict):
            v = _diff_dict(orig.get(k, {}), v)
            if len(v) > 0:
                final[k] = v
        elif v != orig.get(k):
            final[k] = v
    for k, v in orig.items():
        if k not in new:
            final[k] = None
    return final

def _pick_lead_item(items):
    """Pick single representative sample for batch calling to attach calls to.

    For cancer samples, attach to tumor.
    """
    if vcfutils.is_paired_analysis([dd.get_align_bam(x) for x in items], items):
        for data in items:
            if vcfutils.get_paired_phenotype(data) == "tumor":
                return data
        raise ValueError("Did not find tumor sample in paired tumor/normal calling")
    else:
        return items[0]

def get_orig_items(base):
    """Retrieve original items from a diffed set of nested samples.
    """
    assert "group_orig" in base
    out = []
    for data_diff in base["group_orig"]:
        new = utils.deepish_copy(base)
        new.pop("group_orig")
        out.append(_patch_dict(data_diff, new))
    return out

def _patch_dict(diff, base):
    """Patch a dictionary, substituting in changed items from the nested diff.
    """
    for k, v in diff.items():
        if isinstance(v, dict):
            base[k] = _patch_dict(v, base.get(k, {}))
        elif not v:
            base.pop(k, None)
        else:
            base[k] = v
    return base

# ## Split batched variants

def split_variants_by_sample(data):
    """Split a multi-sample call file into inputs for individual samples.

    For tumor/normal paired analyses, do not split the final file and attach
    it to the tumor input.
    """
    # not split, do nothing
    if "group_orig" not in data:
        return [[data]]
    # cancer tumor/normal
    elif (vcfutils.get_paired_phenotype(data)
            and "tumor" in [vcfutils.get_paired_phenotype(d) for d in get_orig_items(data)]):
        out = []
        for i, sub_data in enumerate(get_orig_items(data)):
            if vcfutils.get_paired_phenotype(sub_data) == "tumor":
                cur_batch = tz.get_in(["metadata", "batch"], data)
                if cur_batch:
                    sub_data["metadata"]["batch"] = cur_batch
                sub_data["vrn_file"] = data["vrn_file"]
            else:
                sub_data.pop("vrn_file", None)
            out.append([sub_data])
        return out
    # joint calling or population runs, do not split back up and keep in batches
    else:
        out = []
        for sub_data in get_orig_items(data):
            cur_batch = tz.get_in(["metadata", "batch"], data)
            if cur_batch:
                sub_data["metadata"]["batch"] = cur_batch
            sub_data["vrn_file_batch"] = data["vrn_file"]
            sub_data["vrn_file"] = data["vrn_file"]
            out.append([sub_data])
        return out
