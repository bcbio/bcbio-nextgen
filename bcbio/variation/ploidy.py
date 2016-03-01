"""Calculate expected ploidy for a genomic regions.

Handles configured ploidy, with custom handling for sex chromosomes and pooled
haploid mitochondrial DNA.
"""
import re

import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
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

def _configured_ploidy_sex(items):
    ploidies = set([tz.get_in(["config", "algorithm", "ploidy"], data, 2) for data in items])
    assert len(ploidies) == 1, "Multiple ploidies set for group calling: %s" % ploidies
    ploidy = ploidies.pop()
    sexes = set([str(tz.get_in(["metadata", "sex"], data, "")).lower() for data in items])
    return ploidy, sexes

def get_ploidy(items, region):
    """Retrieve ploidy of a region, handling special cases.
    """
    chrom = chromosome_special_cases(region[0] if isinstance(region, (list, tuple))
                                     else None)
    ploidy, sexes = _configured_ploidy_sex(items)
    if chrom == "mitochondrial":
        # For now, do haploid calling. Could also do pooled calling
        # but not entirely clear what the best default would be.
        return 1
    elif chrom == "X":
        # Do standard diploid calling if we have any females or unspecified.
        if "female" in sexes or "f" in sexes:
            return 2
        elif "male" in sexes or "m" in sexes:
            return 1
        else:
            return 2
    elif chrom == "Y":
        # Always call Y single. If female, filter_vcf_by_sex removes Y regions.
        return 1
    else:
        return ploidy

def filter_vcf_by_sex(vcf_file, items):
    """Post-filter a single sample VCF, handling sex chromosomes.

    Removes Y chromosomes from batches with all female samples.
    """
    out_file = "%s-ploidyfix%s" % utils.splitext_plus(vcf_file)
    if not utils.file_exists(out_file):
        genders = list(_configured_ploidy_sex(items)[-1])
        is_female = len(genders) == 1 and genders[0] and genders[0] in ["female", "f"]
        if is_female:
            orig_out_file = out_file
            out_file = orig_out_file.replace(".vcf.gz", ".vcf")
            with file_transaction(items[0], out_file) as tx_out_file:
                with open(tx_out_file, "w") as out_handle:
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
