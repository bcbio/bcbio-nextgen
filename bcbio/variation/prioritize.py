"""Prioritization scheme for identifying follow up variants in tumor-only samples.

Generalizes the filtering scheme used in VarDict post-processing:

https://github.com/AstraZeneca-NGS/VarDict/blob/9ffec9168e91534fac5fb74b3ec7bdd2badd3464/vcf2txt.pl#L190

The goal is to build up a standard set of prioritization filters based on known
data. Uses GEMINI to load a database of variants with associated third party
query information. Makes use of ExAC, dbSNP, 1000 genomes, clinvar, cosmic and
effects annotations. The general idea is to prioritize deleterious variants
missing or present at a low frequency in the population, or secondarily identified
in external databases like COSMIC and ClinVar.
"""
import csv

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import population, vcfutils
from bcbio.variation import multi as vmulti

def handle_vcf_calls(vcf_file, data):
    """Prioritize VCF calls based on external annotations supplied through GEMINI.
    """
    if not _do_prioritize(data):
        return vcf_file
    else:
        if population.do_db_build([data]):
            gemini_db = population.create_gemini_db(vcf_file, data)
            if gemini_db:
                priority_file = _prep_priority_filter(gemini_db, data)
                return _apply_priority_filter(vcf_file, priority_file, data)
        # No GEMINI database for filtering, return original file
        return vcf_file

def _apply_priority_filter(in_file, priority_file, data):
    """Annotate variants with priority information and use to apply filters.
    """
    out_file = "%s-priority%s" % utils.splitext_plus(in_file)
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            header = ('##INFO=<ID=EPR,Number=.,Type=String,'
                      'Description="Prioritization based on external annotations">')
            header_file = "%s-repeatheader.txt" % utils.splitext_plus(tx_out_file)[0]
            with open(header_file, "w") as out_handle:
                out_handle.write(header)
            cmd = ("bcftools annotate -a {priority_file} -h {header_file} "
                   "-c CHROM,FROM,TO,REF,ALT,INFO/EPR {in_file} | "
                   "bcftools filter -m '+' -s 'LowPriority' "
                   "-e 'EPR[*] != \"pass\"' | "
                   r"""sed 's/\\\"pass\\\"/pass/' | bgzip -c > {tx_out_file}""")
            do.run(cmd.format(**locals()), "Run external annotation based prioritization filtering")
    vcfutils.bgzip_and_index(out_file, data["config"])
    return out_file

def _prep_priority_filter(gemini_db, data):
    """Prepare tabix indexed file with priority based filters and supporting information
    """
    from gemini import GeminiQuery
    out_file = "%s-priority.tsv" % utils.splitext_plus(gemini_db)[0]
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            gq = GeminiQuery(gemini_db)
            attrs = ("chrom, start, end, ref, alt, impact_so, impact_severity, in_dbsnp, "
                     "aaf_esp_all, aaf_1kg_all, aaf_adj_exac_all, cosmic_ids, "
                     "clinvar_sig, clinvar_origin, gt_ref_depths, gt_alt_depths").split(", ")
            gq.run("SELECT %s FROM variants" % ", ".join(attrs))
            sidx = gq.sample_to_idx[dd.get_sample_name(data)]
            header = attrs[:5] + ["filter"] + attrs[5:-2] + ["freq"]
            with open(tx_out_file, "w") as out_handle:
                writer = csv.writer(out_handle, dialect="excel-tab")
                cheader = header[:]
                cheader[0] = "#" + cheader[0]
                writer.writerow(cheader)
                for row in gq:
                    ref_depth = row.gt_ref_depths[sidx]
                    alt_depth = row.gt_alt_depths[sidx]
                    row.row["freq"] = "%.2f" % (float(alt_depth) / float(ref_depth + alt_depth))
                    row.row["filter"] = _calc_priority_filter(row)
                    out = [row[x] for x in header]
                    writer.writerow(out)
    return vcfutils.bgzip_and_index(out_file, data["config"],
                                    tabix_args="-0 -c '#' -s 1 -b 2 -e 3")

def _calc_priority_filter(row):
    """Calculate the priority filter based on external associated data.
    """
    filters = []
    passes = []
    if row["impact_severity"] in ["LOW"]:
        filters.append("lowseverity")
    passes.extend(_find_known(row))
    filters.extend(_known_populations(row))
    if len(filters) == 0 or len(passes) > 0:
        passes.insert(0, "pass")
    return ",".join(passes + filters)

def _known_populations(row):
    """Find variants present in higher frequency in population databases.
    """
    cutoff = 0.1
    out = []
    for pop, key in [("esp", "aaf_esp_all"), ("1000g", "aaf_1kg_all"),
                     ("exac", "aaf_adj_exac_all")]:
        val = row[key]
        if val and val > cutoff:
            out.append(pop)
    return out

def _find_known(row):
    """Find variant present in known pathogenic databases.
    """
    out = []
    clinvar_no = set(["unknown", "untested", "non-pathogenic", "probable-non-pathogenic"])
    if row["cosmic_ids"]:
        out.append("cosmic")
    if (row["clinvar_sig"] and not row["clinvar_sig"] in clinvar_no):
        out.append("clinvar")
    return out

def _do_prioritize(data):
    """Determine if we should perform prioritization.

    Currently done on tumor-only input samples.
    """
    if vcfutils.get_paired_phenotype(data):
        has_tumor = False
        has_normal = False
        for i, sub_data in enumerate(vmulti.get_orig_items(data)):
            if vcfutils.get_paired_phenotype(sub_data) == "tumor":
                has_tumor = True
            elif vcfutils.get_paired_phenotype(sub_data) == "normal":
                has_normal = True
        return has_tumor and not has_normal
