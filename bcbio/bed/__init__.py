import pybedtools as bt
import six

def concat(bed_files, catted=None):
    """
    recursively concat a set of BED files, returning a
    sorted bedtools object of the result
    """
    bed_files = [x for x in bed_files if x]
    if len(bed_files) == 0:
        if catted:
            # move to a .bed extension for downstream tools if not already
            sorted_bed = catted.sort()
            if not sorted_bed.fn.endswith(".bed"):
                return sorted_bed.moveto(sorted_bed.fn + ".bed")
            else:
                return sorted_bed
        else:
            return catted

    if not catted:
        bed_files = list(bed_files)
        catted = bt.BedTool(bed_files.pop())
    else:
        catted = catted.cat(bed_files.pop(), postmerge=False,
                            force_truncate=False)

    return concat(bed_files, catted)

def merge(bedfiles):
    """
    given a BED file or list of BED files merge them an return a bedtools object
    """
    if isinstance(bedfiles, list):
        catted = concat(bedfiles)
    else:
        catted = concat([bedfiles])
    if catted:
        return concat(bedfiles).sort().merge()
    else:
        return catted

def minimize(bed_file):
    """
    strip a BED file down to its three necessary columns: chrom start end
    """
    if not bed_file:
        return bed_file
    else:
        sorted_bed = bt.BedTool(bed_file).cut(range(3)).sort()
        if not sorted_bed.fn.endswith(".bed"):
            return sorted_bed.moveto(sorted_bed.fn + ".bed")
        else:
            return sorted_bed
