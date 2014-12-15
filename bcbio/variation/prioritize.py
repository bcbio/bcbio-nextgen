"""Prioritization scheme for identifying follow up variants in tumor-only samples.

Generalizes the filtering scheme used in VarDict post-processing:

https://github.com/AstraZeneca-NGS/VarDict/blob/9ffec9168e91534fac5fb74b3ec7bdd2badd3464/vcf2txt.pl#L190

The goal is to build up a standard set of prioritization filters based on known
data. Uses GEMINI to load a database of variants with associated third party
query information. Makes use of ExAC, dbSNP, 1000 genomes, clinvar, cosmic and
effects annotations. The general idea is to prioritize deleterious variants
missing or present at a low frequency in the population, or secondarily identified
in external databases.
"""
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
                # TODO filter with information from GEMINI database
                pass
        return vcf_file

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
