"""Extract germline calls from tumor/normal pairs into separate VCF file.

In tumor/normal pipelines, pre-existing germline calls are often of interest
in addition to somatic variants. Different callers distinguish germline calls
in different ways. This unifies the output and extracts into a separate VCF
with germline calls included.
"""
import contextlib

import cyvcf2

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.variation import vcfutils

def extract(data, items):
    """Extract germline calls for the given sample, if tumor/normal or prioritized.
    """
    if vcfutils.get_paired_phenotype(data):
        is_paired = dd.get_batches(data) and len(items) > 1
        if is_paired:
            germline_vcf = _extract_germline(data["vrn_file"], data)
        else:
            germline_vcf = _remove_prioritization(data["vrn_file"], data)
        germline_vcf = vcfutils.bgzip_and_index(germline_vcf, data["config"])
        data["vrn_file_plus"] = {"germline": germline_vcf}
    return data

def _remove_prioritization(in_file, data):
    """Remove tumor-only prioritization and return non-filtered calls.
    """
    out_file = "%s-germline.vcf" % utils.splitext_plus(in_file)[0]
    if not utils.file_uptodate(out_file, in_file) and not utils.file_uptodate(out_file + ".gz", in_file):
        with file_transaction(data, out_file) as tx_out_file:
            reader = cyvcf2.VCF(in_file)
            reader.add_filter_to_header({'ID': 'Somatic', 'Description': 'Variant called as Somatic'})
            # with contextlib.closing(cyvcf2.Writer(tx_out_file, reader)) as writer:
            with open(tx_out_file, "w") as out_handle:
                out_handle.write(reader.raw_header)
                for rec in reader:
                    rec = _update_prioritization_filters(rec)
                    out_handle.write(str(rec))
                    # writer.write_record(rec)
    return out_file

def _update_prioritization_filters(rec):
    rec = _remove_filter(rec, "LowPriority")
    return _update_germline_filters(rec)

def _extract_germline(in_file, data):
    """Extract germline calls non-somatic, non-filtered calls.
    """
    out_file = "%s-germline.vcf" % utils.splitext_plus(in_file)[0]
    if not utils.file_uptodate(out_file, in_file) and not utils.file_uptodate(out_file + ".gz", in_file):
        with file_transaction(data, out_file) as tx_out_file:
            reader = cyvcf2.VCF(in_file)
            reader.add_filter_to_header({'ID': 'Somatic', 'Description': 'Variant called as Somatic'})
            #with contextlib.closing(cyvcf2.Writer(tx_out_file, reader)) as writer:
            with open(tx_out_file, "w") as out_handle:
                out_handle.write(reader.raw_header)
                for rec in reader:
                    rec = _update_germline_filters(rec)
                    out_handle.write(str(rec))
                    #writer.write_record(rec)
    return out_file

def _update_germline_filters(rec):
    rec = _remove_germline_filter(rec, "REJECT")
    rec = _remove_germline_filter(rec, "germline_risk")
    rec = _add_somatic_filter(rec)
    return rec

def _add_somatic_filter(rec):
    if _is_somatic(rec):
        return _add_filter(rec, "Somatic")
    return rec

def _remove_germline_filter(rec, name):
    """Check if germline based on STATUS/SS and REJECT flag.

    Handles VarDict, FreeBayes, MuTect, MuTect2 and VarScan.
    """
    if _is_germline(rec):
        if rec.FILTER and name in rec.FILTER:
            return _remove_filter(rec, name)
    elif not _is_somatic(rec):
        if rec.FILTER and name in rec.FILTER:
            return _remove_filter(rec, name)
    return rec

def _is_somatic(rec):
    """Handle somatic classifications from MuTect, MuTect2, VarDict and VarScan
    """
    if _has_somatic_flag(rec):
        return True
    if _is_mutect2_somatic(rec):
        return True
    ss_flag = rec.INFO.get("SS")
    if ss_flag is not None:
        if str(ss_flag) == "2":
            return True
    status_flag = rec.INFO.get("STATUS")
    if status_flag is not None:
        if str(status_flag).lower() in ["somatic", "likelysomatic", "strongsomatic", "samplespecific"]:
            return True
    return False

def _has_somatic_flag(rec):
    try:
        rec.INFO["SOMATIC"]
        return True
    except KeyError:
        return False

def _is_mutect2_somatic(rec):
    """MuTect2 does not use SOMATIC flag, instead using presence of tumor TLOD and PASS.
    """
    return rec.INFO.get("TLOD") is not None and rec.FILTER is None

def _is_germline(rec):
    """Handle somatic INFO classifications from MuTect, MuTect2, VarDict and VarScan
    """
    if _has_somatic_flag(rec):
        return False
    if _is_mutect2_somatic(rec):
        return False
    ss_flag = rec.INFO.get("SS")
    if ss_flag is not None:
        if str(ss_flag) == "1":
            return True
    status_flag = rec.INFO.get("STATUS")
    if status_flag is not None:
        if str(status_flag).lower() in ["germline", "likelyloh", "strongloh", "afdiff", "deletion"]:
            return True
    return False

def _add_filter(rec, name):
    if rec.FILTER:
        filters = rec.FILTER.split(";")
    else:
        filters = []
    if name not in filters:
        filters.append(name)
        rec.FILTER = filters
    return rec

def _remove_filter(rec, name):
    """Remove filter with the given name from VCF input.
    """
    if rec.FILTER:
        filters = rec.FILTER.split(";")
    else:
        filters = []
    new_filters = [x for x in filters if not str(x) == name]
    if len(new_filters) == 0:
        new_filters = ["PASS"]
    rec.FILTER = new_filters
    return rec
