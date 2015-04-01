import pybedtools as bt

def concat(bed_files, catted=None):
    """
    recursively concat a set of BED files, returning a
    sorted bedtools object of the result
    """
    if len(bed_files) == 0:
        return catted.sort()

    if not catted:
        bed_files = list(bed_files)
        catted = bt.BedTool(bed_files.pop())
    else:
        catted = catted.cat(bed_files.pop(), postmerge=False,
                            force_truncate=False)

    return concat(bed_files, catted)
