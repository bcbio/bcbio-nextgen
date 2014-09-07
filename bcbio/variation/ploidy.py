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
    sexes = set([tz.get_in(["metadata", "sex"], data, "").lower() for data in items])
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
        if "female" in sexes:
            return 2
        elif "male" in sexes:
            return 1
        else:
            return 2
    elif chrom == "Y":
        # Always call Y single. If female, filter_vcf_by_sex removes Y regions.
        return 1
    else:
        return ploidy

def _to_haploid(parts):
    """Check if a variant call is homozygous variant, convert to haploid.
    XXX Needs generalization or use of a standard VCF library.
    """
    finfo = dict(zip(parts[-2].split(":"), parts[-1].strip().split(":")))
    pat = re.compile(r"\||/")
    if "GT" in finfo:
        calls = set(pat.split(finfo["GT"]))
        if len(calls) == 1:
            gt_index = parts[-2].split(":").index("GT")
            call_parts = parts[-1].strip().split(":")
            call_parts[gt_index] = calls.pop()
            parts[-1] = ":".join(call_parts) + "\n"
            return "\t".join(parts)

def _fix_line_ploidy(line, sex):
    """Check variant calls to be sure if conforms to expected ploidy for sex/custom chromosomes.
    """
    parts = line.split("\t")
    chrom = chromosome_special_cases(parts[0])
    if chrom == "mitochondrial":
        return _to_haploid(parts)
    elif chrom == "X":
        if sex == "male":
            return _to_haploid(parts)
        else:
            return line
    elif chrom == "Y":
        if sex != "female":
            return _to_haploid(parts)
    else:
        return line

def filter_vcf_by_sex(vcf_file, data):
    """Post-filter a single sample VCF, handling sex chromosomes.

    Handles sex chromosomes and mitochondrial. Does not try to resolve called
    hets into potential homozygotes when converting diploid to haploid.

    Skips filtering on pooled samples, we still need to implement.
    """
    if len(vcfutils.get_samples(vcf_file)) > 1:
        return vcf_file
    _, sexes = _configured_ploidy_sex([data])
    sex = sexes.pop()
    out_file = "%s-ploidyfix%s" % utils.splitext_plus(vcf_file)
    if not utils.file_exists(out_file):
        orig_out_file = out_file
        out_file = orig_out_file.replace(".vcf.gz", ".vcf")
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                with utils.open_gzipsafe(vcf_file) as in_handle:
                    for line in in_handle:
                        if line.startswith("#"):
                            out_handle.write(line)
                        else:
                            line = _fix_line_ploidy(line, sex)
                            if line:
                                out_handle.write(line)
        if orig_out_file.endswith(".gz"):
            out_file = vcfutils.bgzip_and_index(out_file, data["config"])
    return out_file
