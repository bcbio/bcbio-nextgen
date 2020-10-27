"""Calculate expected ploidy for a genomic regions.

Handles configured ploidy, with custom handling for sex chromosomes and pooled
haploid mitochondrial DNA.
"""
import collections
import io

import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.variation import vcfutils

def chromosome_special_cases(chrom):
    if chrom in ["MT", "M", "chrM", "chrMT"]:
        return "mitochondrial"
    elif chrom in ["X", "chrX"]:
        return "X"
    elif chrom in ["Y", "chrY"]:
        return "Y"
    else:
        return chrom

def _configured_ploidy(items):
    ploidies = collections.defaultdict(set)
    for data in items:
        ploidy = dd.get_ploidy(data)
        if isinstance(ploidy, dict):
            for k, v in ploidy.items():
                ploidies[k].add(v)
        else:
            ploidies["default"].add(ploidy)
    out = {}
    for k, vs in ploidies.items():
        assert len(vs) == 1, "Multiple ploidies set for group calling: %s %s" % (k, list(vs))
        out[k] = vs.pop()
    return out

def _configured_genders(items):
    return set([str(tz.get_in(["metadata", "sex"], data, "")).lower() for data in items])

def get_ploidy(items, region=None):
    """Retrieve ploidy of a region, handling special cases.
    """
    chrom = chromosome_special_cases(region[0] if isinstance(region, (list, tuple))
                                     else None)
    ploidy = _configured_ploidy(items)
    sexes = _configured_genders(items)
    if chrom == "mitochondrial":
        # For now, do haploid calling. Could also do pooled calling
        # but not entirely clear what the best default would be.
        return ploidy.get("mitochondrial", 1)
    elif chrom == "X":
        # Do standard diploid calling if we have any females or unspecified.
        if "female" in sexes or "f" in sexes:
            return ploidy.get("female", ploidy["default"])
        elif "male" in sexes or "m" in sexes:
            return ploidy.get("male", 1)
        else:
            return ploidy.get("female", ploidy["default"])
    elif chrom == "Y":
        # Always call Y single. If female, filter_vcf_by_sex removes Y regions.
        return 1
    else:
        return ploidy["default"]

def filter_vcf_by_sex(vcf_file, items):
    """Post-filter a single sample VCF, handling sex chromosomes.

    Removes Y chromosomes from batches with all female samples.
    """
    out_file = "%s-ploidyfix%s" % utils.splitext_plus(vcf_file)
    if not utils.file_exists(out_file):
        genders = list(_configured_genders(items))
        is_female = len(genders) == 1 and genders[0] and genders[0] in ["female", "f"]
        if is_female:
            orig_out_file = out_file
            out_file = orig_out_file.replace(".vcf.gz", ".vcf")
            with file_transaction(items[0], out_file) as tx_out_file:
                with io.open(tx_out_file, "w", encoding="utf-8") as out_handle:
                    with utils.open_gzipsafe(vcf_file) as in_handle:
                        for line in in_handle:
                            if line.startswith("#"):
                                out_handle.write(line)
                            else:
                                chrom = chromosome_special_cases(line.split("\t"))
                                if chrom != "Y":
                                    out_handle.write(line)
            if orig_out_file.endswith(".gz"):
                out_file = vcfutils.bgzip_and_index(out_file, items[0]["config"])
        else:
            out_file = vcf_file
    return out_file
