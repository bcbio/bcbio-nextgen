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
from bcbio import bed

def breakpoints_by_caller(bed_files):
    """
    given a list of BED files of the form
    chrom start end caller
    return a BedTool of breakpoints as each line with
    the fourth column the caller with evidence for the breakpoint
    chr1 1 10 caller1 -> chr1 1 1 caller1
    chr1 1 20 caller2    chr1 1 1 caller2
                         chr1 10 10 caller1
                         chr1 20 20 caller2
    """
    merged = concat(bed_files)
    if not merged:
        return []
    grouped_start = merged.groupby(g=[1, 2, 2], c=4, o=["distinct"]).filter(lambda r: r.end > r.start).saveas()
    grouped_end = merged.groupby(g=[1, 3, 3], c=4, o=["distinct"]).filter(lambda r: r.end > r.start).saveas()
    together = concat([grouped_start, grouped_end])
    if together:
        final = together.expand(c=4)
        final = final.sort()
        return final

def _get_sv_callers(items):
    """
    return a sorted list of all of the structural variant callers run
    """
    callers = []
    for data in items:
        for sv in data.get("sv", []):
            callers.append(sv["variantcaller"])
    return list(set([x for x in callers if x != "sv-ensemble"])).sort()

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
    out = []
    added = False
    for data in items:
        for sv in data.get("sv", []):
            if not added and sv["variantcaller"] == "sv-ensemble":
                added = True
                if "plot" not in sv:
                    sv["plot"] = {}
                sv["plot"]["coverage"] = plot
        out.append(data)
    return out

def _prioritize_plot_regions(region_bt, data):
    """Avoid plotting large numbers of regions due to speed issues. Prioritize most interesting.

    XXX For now, just removes larger regions and avoid plotting thousands of regions.
    Longer term we'll insert biology-based prioritization.
    """
    max_plots = 1000
    max_size = 100 * 1000 # 100kb
    out_file = "%s-priority%s" % utils.splitext_plus(region_bt.fn)
    num_plots = 0
    if not utils.file_uptodate(out_file, region_bt.fn):
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                for r in region_bt:
                    if r.stop - r.start < max_size:
                        if num_plots < max_plots:
                            num_plots += 1
                            out_handle.write(str(r))
    return out_file

def by_regions(items):
    """Plot for a union set of combined ensemble regions across all of the data
       items.
    """
    work_dir = os.path.join(dd.get_work_dir(items[0]), "structural", "coverage")
    safe_makedir(work_dir)
    out_file = os.path.join(work_dir, "%s-coverage.pdf" % (dd.get_sample_name(items[0])))
    if file_exists(out_file):
        items = _add_regional_coverage_plot(items, out_file)
    else:
        bed_files = _get_ensemble_bed_files(items)
        merged = bed.merge(bed_files)
        breakpoints = breakpoints_by_caller(bed_files)
        if merged:
            priority_merged = _prioritize_plot_regions(merged, items[0])
            out_file = plot_multiple_regions_coverage(items, out_file,
                                                      priority_merged, breakpoints)
            items = _add_regional_coverage_plot(items, out_file)
    return items
