"""Provide plots of structural variations to manually validate results.
Uses existing plots from CNVkit along with custom plotting of coverage to
provide the ability to quickly validate and explore predicted structural
variants.
"""
import os
from bcbio.pipeline import datadict as dd
from bcbio.variation import vcfutils
from bcbio.bam.coverage import plot_multiple_regions_coverage
from bcbio.bed import concat
from bcbio.utils import file_exists, safe_makedir

def _merge_sv_calls(bed_files):
    """
    merge a set of structural variant BED files and return a bedtools object
    """
    if bed_files:
        merged = concat(bed_files)
        merged = merged.merge()
        return(merged)

def _get_ensemble_bed_files(items):
    """
    get all ensemble structural BED file calls, skipping any normal samples from
    tumor/normal calls
    """
    bed_files = []
    for data in items:
        for sv in data.get("sv", []):
            if sv["variantcaller"] == "sv-ensemble":
                if ("vrn_file" in sv and not vcfutils.get_paired_phenotype(data) == "normal"
                      and file_exists(sv["vrn_file"])):
                    bed_files.append(sv["vrn_file"])
    return bed_files

def _add_regional_coverage_plot(items, plot):
    for data in items:
        for sv in data.get("sv", []):
            if sv["variantcaller"] == "sv-ensemble":
                sv["coverage_plot"] = plot
    return items

def by_regions(items):
    """Plot for a union set of combined ensemble regions across all of the data
       items.
    """
    work_dir = os.path.join(dd.get_work_dir(items[0]), "structural", "coverage")
    safe_makedir(work_dir)
    out_file = os.path.join(work_dir, "coverage.pdf")
    if file_exists(out_file):
        items = _add_regional_coverage_plot(items, out_file)
    else:
        merged = _merge_sv_calls(_get_ensemble_bed_files(items))
        if merged:
            out_file = plot_multiple_regions_coverage(items, out_file, merged)
            items = _add_regional_coverage_plot(items, out_file)
    return items
