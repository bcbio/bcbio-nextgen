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
import collections
import csv
import re

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import population, vcfutils

geneimpacts = utils.LazyImport("geneimpacts")
cyvcf2 = utils.LazyImport("cyvcf2")

def handle_vcf_calls(vcf_file, data, orig_items):
    """Prioritize VCF calls based on external annotations supplied through GEMINI.
    """
    if not _do_prioritize(orig_items):
        return vcf_file
    else:
        ann_vcf = population.run_vcfanno(vcf_file, data)
        if ann_vcf:
            priority_file = _prep_priority_filter_vcfanno(ann_vcf, data)
            return _apply_priority_filter(ann_vcf, priority_file, data)
        # No data available for filtering, return original file
        else:
            return vcf_file

def _apply_priority_filter(in_file, priority_file, data):
    """Annotate variants with priority information and use to apply filters.
    """
    out_file = "%s-priority%s" % utils.splitext_plus(in_file)
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            header = ('##INFO=<ID=EPR,Number=.,Type=String,'
                      'Description="Somatic prioritization based on external annotations, '
                      'identify as likely germline">')
            header_file = "%s-repeatheader.txt" % utils.splitext_plus(tx_out_file)[0]
            with open(header_file, "w") as out_handle:
                out_handle.write(header)
            if "tumoronly_germline_filter" in dd.get_tools_on(data):
                filter_cmd = ("bcftools filter -m '+' -s 'LowPriority' "
                              """-e "EPR[0] != 'pass'" |""")
            else:
                filter_cmd = ""
            cmd = ("bcftools annotate -a {priority_file} -h {header_file} "
                   "-c CHROM,FROM,TO,REF,ALT,INFO/EPR {in_file} | "
                   "{filter_cmd} bgzip -c > {tx_out_file}")
            do.run(cmd.format(**locals()), "Run external annotation based prioritization filtering")
    vcfutils.bgzip_and_index(out_file, data["config"])
    return out_file

def _prep_priority_filter_vcfanno(in_vcf, data):
    """Prepare tabix file with priority filters based on vcfanno annotations.
    """
    pops = ['af_adj_exac_afr', 'af_adj_exac_amr', 'af_adj_exac_eas',
            'af_adj_exac_fin', 'af_adj_exac_nfe', 'af_adj_exac_oth', 'af_adj_exac_sas',
            'af_exac_all', 'max_aaf_all',
            "af_esp_ea", "af_esp_aa", "af_esp_all", "af_1kg_amr", "af_1kg_eas",
            "af_1kg_sas", "af_1kg_afr", "af_1kg_eur", "af_1kg_all"]
    known = ["cosmic_ids", "cosmic_id", "clinvar_sig"]
    out_file = "%s-priority.tsv" % utils.splitext_plus(in_vcf)[0]
    if not utils.file_exists(out_file) and not utils.file_exists(out_file + ".gz"):
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                writer = csv.writer(out_handle, dialect="excel-tab")
                header = ["#chrom", "start", "end", "ref", "alt", "filter"]
                writer.writerow(header)
                vcf_reader = cyvcf2.VCF(in_vcf)
                impact_info = _get_impact_info(vcf_reader)
                for rec in vcf_reader:
                    row = _prepare_vcf_rec(rec, pops, known, impact_info)
                    cur_filter = _calc_priority_filter(row, pops)
                    writer.writerow([rec.CHROM, rec.start, rec.end, rec.REF, ",".join(rec.ALT), cur_filter])
    return vcfutils.bgzip_and_index(out_file, data["config"],
                                    tabix_args="-0 -c '#' -s 1 -b 2 -e 3")

def _get_impact_info(vcf_reader):
    """Retrieve impact parsing information from INFO header.
    """
    ImpactInfo = collections.namedtuple("ImpactInfo", "header, gclass, id")
    KEY_2_CLASS = {
        'CSQ': geneimpacts.VEP,
        'ANN': geneimpacts.SnpEff,
        'BCSQ': geneimpacts.BCFT}
    for l in (x.strip() for x in _from_bytes(vcf_reader.raw_header).split("\n")):
        if l.startswith("##INFO"):
            patt = re.compile(r"(\w+)=(\"[^\"]+\"|[^,]+)")
            stub = l.split("=<")[1].rstrip(">")
            d = dict(patt.findall(_from_bytes(stub)))
            if d["ID"] in KEY_2_CLASS:
                return ImpactInfo(_parse_impact_header(d), KEY_2_CLASS[d["ID"]], d["ID"])

def _from_bytes(s):
    if isinstance(s, bytes):
        import locale
        ENC = locale.getpreferredencoding()
        try:
            return s.decode(ENC)
        except UnicodeDecodeError:
            return s.decode('utf8')
    return s

def _parse_impact_header(hdr_dict):
    """Parse fields for impact, taken from vcf2db
    """
    desc = hdr_dict["Description"]
    if hdr_dict["ID"] == "ANN":
        parts = [x.strip("\"'") for x in re.split("\s*\|\s*", desc.split(":", 1)[1].strip('" '))]
    elif hdr_dict["ID"] == "EFF":
        parts = [x.strip(" [])'(\"") for x in re.split("\||\(", desc.split(":", 1)[1].strip())]
    elif hdr_dict["ID"] == "CSQ":
        parts = [x.strip(" [])'(\"") for x in re.split("\||\(", desc.split(":", 1)[1].strip())]
    elif hdr_dict["ID"] == "BCSQ":
        parts = desc.split(']', 1)[1].split(']')[0].replace('[','').split("|")
    else:
        raise Exception("don't know how to use %s as annotation" % hdr_dict["ID"])
    return parts

def _prepare_vcf_rec(rec, pops, known, impact_info):
    """Parse a vcfanno output into a dictionary of useful attributes.
    """
    out = {}
    for k in pops + known:
        out[k] = rec.INFO.get(k)
    if impact_info:
        cur_info = rec.INFO.get(impact_info.id)
        if cur_info:
            cur_impacts = [impact_info.gclass(e, impact_info.header) for e in _from_bytes(cur_info).split(",")]
            top = geneimpacts.Effect.top_severity(cur_impacts)
            if isinstance(top, list):
                top = top[0]
            out["impact_severity"] = top.effect_severity
    return out

def _calc_priority_filter(row, pops):
    """Calculate the priority filter based on external associated data.

    - Pass high/medium impact variants not found in population databases
    - Pass variants found in COSMIC or Clinvar provided they don't have two
      additional reasons to filter (found in multiple external populations)
    """
    filters = []
    passes = []
    passes.extend(_find_known(row))
    filters.extend(_known_populations(row, pops))
    if len(filters) == 0 or (len(passes) > 0 and len(filters) < 2):
        passes.insert(0, "pass")
    return ",".join(passes + filters)

def _known_populations(row, pops):
    """Find variants present in substantial frequency in population databases.
    """
    cutoff = 0.01
    out = set([])
    for pop, base in [("esp", "af_esp_all"), ("1000g", "af_1kg_all"),
                      ("exac", "af_exac_all"), ("anypop", "max_aaf_all")]:
        for key in [x for x in pops if x.startswith(base)]:
            val = row[key]
            if val and val > cutoff:
                out.add(pop)
    return sorted(list(out))

def _find_known(row):
    """Find variant present in known pathogenic databases.
    """
    out = []
    clinvar_no = set(["unknown", "untested", "non-pathogenic", "probable-non-pathogenic",
                      "uncertain_significance", "uncertain_significance", "not_provided",
                      "benign", "likely_benign"])
    if row["cosmic_ids"] or row["cosmic_id"]:
        out.append("cosmic")
    if row["clinvar_sig"] and not row["clinvar_sig"].lower() in clinvar_no:
        out.append("clinvar")
    return out

def _do_prioritize(items):
    """Determine if we should perform prioritization.

    Currently done on tumor-only input samples and feeding into PureCN
    which needs the germline annotations.
    """
    if not any("tumoronly-prioritization" in dd.get_tools_off(d) for d in items):
        if vcfutils.get_paired_phenotype(items[0]):
            has_tumor = False
            has_normal = False
            for sub_data in items:
                if vcfutils.get_paired_phenotype(sub_data) == "tumor":
                    has_tumor = True
                elif vcfutils.get_paired_phenotype(sub_data) == "normal":
                    has_normal = True
            return has_tumor and not has_normal
