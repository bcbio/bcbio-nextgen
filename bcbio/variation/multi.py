"""Organize samples for coordinated multi-sample processing.

Handles grouping of related families or batches to go through variant
calling simultaneously.
"""
import collections
import copy
import os

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction

def group_families(xs):
    """Group samples into families for simultaneous variant calling.

    Identify all samples to call together: those in the same family,
    variant caller and genomic region.
    Pull together all BAM files from this family and process together,
    Provide details to pull these finalized files back into individual
    expected files.
    """
    singles = []
    family_groups = collections.defaultdict(list)
    for data, region, out_fname in xs:
        family = data.get("metadata", {}).get("family")
        caller = data["config"]["algorithm"]["variantcaller"]
        if family is not None:
            family_groups[(family, region, caller)].append((data, out_fname))
        else:
            singles.append((data, region, out_fname))
    families = []
    remap_families = {}
    for (family, region, _), xs in family_groups.iteritems():
        cur_data, cur_fname = xs[0]
        family_fname = utils.append_stem(cur_fname, family, "-")
        family_data = copy.deepcopy(cur_data)
        family_data["work_bam"] = [x[0]["work_bam"] for x in xs]
        family_data["group"] = family_fname
        families.append((family_data, region, family_fname))
        remap_families[family_fname] = xs
    return singles + families, remap_families

def split_variants_by_sample(data):
    """Split a multi-sample call file into inputs for individual samples.
    """
    config = data["config"]
    vrn_file = data["vrn_file"]
    out = []
    for sub_data, sub_vrn_file in data["group_orig"]:
        if is_multisample(vrn_file):
            select_sample_from_vcf(vrn_file, sub_data["name"][-1], sub_vrn_file,
                                   data["sam_ref"], config)
        else:
            os.symlink(vrn_file, sub_vrn_file)
        sub_data["vrn_file"] = sub_vrn_file
        out.append(sub_data)
    return out

def select_sample_from_vcf(in_file, sample, out_file, ref_file, config):
    """Select a single sample from the supplied multisample VCF file.
    """
    brunner = broad.runner_from_config(config)
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "SelectVariants",
                      "-R", ref_file,
                      "--sample_name", sample,
                      "--variant", in_file,
                      "--out", tx_out_file]
            brunner.run_gatk(params)
    return out_file

def is_multisample(fname):
    """Check VCF header to determine if we have a multi-sample file.
    """
    with open(fname) as in_handle:
        for line in in_handle:
            if line.startswith("#CHROM"):
                return len(line.split("\t")) > 10
