import pybedtools as bt
import six

def concat(bed_files, catted=None):
    """
    recursively concat a set of BED files, returning a
    sorted bedtools object of the result
    """
    if len(bed_files) == 0:
        if catted:
            return catted.sort()
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
