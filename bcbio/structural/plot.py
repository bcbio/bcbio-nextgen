"""Provide plots of structural variations to manually validate results.
Uses existing plots from CNVkit along with custom plotting of coverage to
provide the ability to quickly validate and explore predicted structural
variants.
"""
import os

import pybedtools

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.variation import vcfutils
from bcbio.bam.coverage import plot_multiple_regions_coverage
from bcbio.bed import concat
from bcbio.utils import file_exists, safe_makedir

def _merge_sv_calls(bed_files, out_file, data):
    """
    merge a set of structural variant BED files and return a bedtools object
    """
    if bed_files:
        if not utils.file_exists(out_file):
            with file_transaction(data, out_file) as tx_out_file:
                merged = concat(bed_files)
                merged = merged.sort().merge().saveas(tx_out_file)
        return pybedtools.BedTool(out_file)

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

def _prioritize_plot_regions(region_bt, data):
    """Avoid plotting large numbers of regions due to speed issues. Prioritize most interesting.

    XXX For now, just removes larger regions. Longer term we'll insert biology-based
    prioritization.
    """
    max_size = 100 * 1000 # 100kb
    out_file = "%s-priority%s" % utils.splitext_plus(region_bt.fn)
    if not utils.file_uptodate(out_file, region_bt.fn):
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                for r in region_bt:
                    if r.stop - r.start < max_size:
                        out_handle.write(str(r))
    return out_file

def by_regions(items):
    """Plot for a union set of combined ensemble regions across all of the data
       items.
    """
    work_dir = os.path.join(dd.get_work_dir(items[0]), "structural", "coverage")
    safe_makedir(work_dir)
    out_file = os.path.join(work_dir, "%s-coverage.pdf" % (dd.get_sample_name(items[0])))
    merged_file = "%s-inputs.bed" % utils.splitext_plus(out_file)[0]
    if file_exists(out_file):
        items = _add_regional_coverage_plot(items, out_file)
    else:
        merged = _merge_sv_calls(_get_ensemble_bed_files(items), merged_file, items[0])
        if merged:
            priority_merged = _prioritize_plot_regions(merged, items[0])
            out_file = plot_multiple_regions_coverage(items, out_file, priority_merged)
            items = _add_regional_coverage_plot(items, out_file)
    return items
