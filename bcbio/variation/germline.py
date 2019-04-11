"""Extract germline calls from tumor/normal pairs into separate VCF file.

In tumor/normal pipelines, pre-existing germline calls are often of interest
in addition to somatic variants. Different callers distinguish germline calls
in different ways. This unifies the output and extracts into a separate VCF
with germline calls included.
"""
import os
import contextlib

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import vcfutils

cyvcf2 = utils.LazyImport("cyvcf2")

def split_somatic(items):
    """Split somatic batches, adding a germline target.

    Enables separate germline calling of samples using shared alignments.
    """
    items = [_clean_flat_variantcaller(x) for x in items]
    somatic_groups, somatic, non_somatic = vcfutils.somatic_batches(items)
    # extract germline samples to run from normals in tumor/normal pairs
    germline_added = set([])
    germline = []
    for somatic_group in somatic_groups:
        paired = vcfutils.get_paired(somatic_group)
        if paired and paired.normal_data:
            cur = utils.deepish_copy(paired.normal_data)
            vc = dd.get_variantcaller(cur)
            if isinstance(vc, dict) and "germline" in vc:
                if cur["description"] not in germline_added:
                    germline_added.add(cur["description"])
                    cur["rgnames"]["sample"] = cur["description"]
                    cur["metadata"]["batch"] = "%s-germline" % cur["description"]
                    cur["metadata"]["phenotype"] = "germline"
                    cur = remove_align_qc_tools(cur)
                    cur["config"]["algorithm"]["variantcaller"] = vc["germline"]
                    germline.append(cur)
    # Fix variantcalling specification for only somatic targets
    somatic_out = []
    for data in somatic:
        vc = dd.get_variantcaller(data)
        if isinstance(vc, dict) and "somatic" in vc:
            data["config"]["algorithm"]["variantcaller"] = vc["somatic"]
        somatic_out.append(data)
    return non_somatic + somatic_out + germline

def _clean_flat_variantcaller(data):
    """Convert flattened dictionary from CWL representation into dictionary.

    CWL flattens somatic/germline tags into a set of strings, which we
    reconstitute as a dictionary.
    """
    vc = dd.get_variantcaller(data)
    if isinstance(vc, (list, tuple)) and all([x.count(":") == 1 for x in vc]):
        out = {}
        for v in vc:
            k, v = v.split(":")
            if k in out:
                out[k].append(v)
            else:
                out[k] = [v]
        data = dd.set_variantcaller(data, out)
    return data

def remove_align_qc_tools(data):
    """Remove alignment based QC tools we don't need for data replicates.

    When we do multiple variant calling on a sample file (somatic/germline),
    avoid re-running QC.
    """
    align_qc = set(["qsignature", "coverage", "picard", "samtools", "fastqc"])
    data["config"]["algorithm"]["qc"] = [t for t in dd.get_algorithm_qc(data)
                                         if t not in align_qc]
    return data

def extract(data, items, out_dir=None):
    """Extract germline calls for the given sample, if tumor only.
    """
    if vcfutils.get_paired_phenotype(data):
        if len(items) == 1:
            germline_vcf = _remove_prioritization(data["vrn_file"], data, out_dir)
            germline_vcf = vcfutils.bgzip_and_index(germline_vcf, data["config"])
            data["vrn_file_plus"] = {"germline": germline_vcf}
    return data

def filter_to_pass_and_reject(in_file, paired, out_dir=None):
    """Filter VCF to only those with a strict PASS/REJECT: somatic + germline.

    Removes low quality calls filtered but also labeled with REJECT.
    """
    from bcbio.heterogeneity import bubbletree
    out_file = "%s-prfilter.vcf.gz" % utils.splitext_plus(in_file)[0]
    if out_dir:
        out_file = os.path.join(out_dir, os.path.basename(out_file))
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(paired.tumor_data, out_file) as tx_out_file:
            max_depth = bubbletree.max_normal_germline_depth(in_file, bubbletree.PARAMS, paired)
            tx_out_plain = tx_out_file.replace(".vcf.gz", ".vcf")
            with contextlib.closing(cyvcf2.VCF(in_file)) as reader:
                reader = _add_db_to_header(reader)
                with contextlib.closing(cyvcf2.Writer(tx_out_plain, reader)) as writer:
                    for rec in reader:
                        filters = rec.FILTER.split(";") if rec.FILTER else []
                        other_filters = [x for x in filters if x not in ["PASS", ".", "REJECT"]]
                        if len(other_filters) == 0 or bubbletree.is_info_germline(rec):
                            # Germline, check if we should include based on frequencies
                            if "REJECT" in filters or bubbletree.is_info_germline(rec):
                                stats = bubbletree._is_possible_loh(rec, reader, bubbletree.PARAMS, paired,
                                                                    use_status=True, max_normal_depth=max_depth)
                                if stats:
                                    rec.FILTER = "PASS"
                                    rec.INFO["DB"] = True
                                    writer.write_record(rec)
                            # Somatic, always include
                            else:
                                writer.write_record(rec)
            vcfutils.bgzip_and_index(tx_out_plain, paired.tumor_data["config"])
    return out_file

def _add_db_to_header(reader):
    try:
        reader["DB"]
    except KeyError:
        reader.add_info_to_header({'ID': 'DB', 'Description': 'Likely germline variant',
                                   'Type': 'Flag', 'Number': '0'})
    return reader

def fix_germline_samplename(in_file, sample_name, data):
    """Replace germline sample names, originally from normal BAM file.
    """
    out_file = "%s-fixnames%s" % utils.splitext_plus(in_file)
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            sample_file = "%s-samples.txt" % utils.splitext_plus(tx_out_file)[0]
            with open(sample_file, "w") as out_handle:
                out_handle.write("%s\n" % sample_name)
            cmd = ("bcftools reheader -s {sample_file} {in_file} -o {tx_out_file}")
            do.run(cmd.format(**locals()), "Fix germline samplename: %s" % sample_name)
    return vcfutils.bgzip_and_index(out_file, data["config"])

def _remove_prioritization(in_file, data, out_dir=None):
    """Remove tumor-only prioritization and return non-filtered calls.
    """
    out_file = "%s-germline.vcf" % utils.splitext_plus(in_file)[0]
    if out_dir:
        out_file = os.path.join(out_dir, os.path.basename(out_file))
    if not utils.file_uptodate(out_file, in_file) and not utils.file_uptodate(out_file + ".gz", in_file):
        with file_transaction(data, out_file) as tx_out_file:
            reader = cyvcf2.VCF(str(in_file))
            reader.add_filter_to_header({'ID': 'Somatic', 'Description': 'Variant called as Somatic'})
            # with open(tx_out_file, "w") as out_handle:
            #     out_handle.write(reader.raw_header)
            with contextlib.closing(cyvcf2.Writer(tx_out_file, reader)) as writer:
                for rec in reader:
                    rec = _update_prioritization_filters(rec)
                    # out_handle.write(str(rec))
                    writer.write_record(rec)
    return out_file

def _update_prioritization_filters(rec):
    rec = vcfutils.cyvcf_remove_filter(rec, "LowPriority")
    return _update_germline_filters(rec)

def _extract_germline(in_file, data):
    """Extract germline calls non-somatic, non-filtered calls.
    """
    out_file = "%s-germline.vcf" % utils.splitext_plus(in_file)[0]
    if not utils.file_uptodate(out_file, in_file) and not utils.file_uptodate(out_file + ".gz", in_file):
        with file_transaction(data, out_file) as tx_out_file:
            reader = cyvcf2.VCF(str(in_file))
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
        return vcfutils.cyvcf_add_filter(rec, "Somatic")
    return rec

def _remove_germline_filter(rec, name):
    """Check if germline based on STATUS/SS and REJECT flag.

    Handles VarDict, FreeBayes, MuTect, MuTect2 and VarScan.
    """
    if _is_germline(rec):
        if rec.FILTER and name in rec.FILTER:
            return vcfutils.cyvcf_remove_filter(rec, name)
    elif not _is_somatic(rec):
        if rec.FILTER and name in rec.FILTER:
            return vcfutils.cyvcf_remove_filter(rec, name)
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
    epr = rec.INFO.get("EPR", "").split(",")
    if epr and all([p == "pass" for p in epr]):
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
    """Handle somatic INFO classifications from MuTect, MuTect2, VarDict, VarScan and Octopus.
    """
    if _has_somatic_flag(rec):
        return False
    if _is_mutect2_somatic(rec):
        return False
    ss_flag = rec.INFO.get("SS")
    if ss_flag is not None:
        if str(ss_flag) == "1":
            return True
    # Octopus, assessed for potentially being Germline and not flagged SOMATIC
    # https://github.com/luntergroup/octopus/wiki/Calling-models:-Cancer#qual-vs-pp
    pp = rec.INFO.get("PP")
    if pp and float(pp) / float(rec.QUAL) >= 0.5:
        return True
    status_flag = rec.INFO.get("STATUS")
    if status_flag is not None:
        if str(status_flag).lower() in ["germline", "likelyloh", "strongloh", "afdiff", "deletion"]:
            return True
    return False
