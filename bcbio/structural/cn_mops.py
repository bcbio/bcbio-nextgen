"""Copy number detection using read counts, with cn.mops.

http://www.bioconductor.org/packages/release/bioc/html/cn.mops.html
"""

def run(items):
    """Detect copy number variations from batched set of samples using cn.mops.
    """
    pass

_script = """
bam_files <- strsplit("{bam_file_str}", ",")[[1]]
sample_names <- strsplit("{names_str}", ",")[[1]]

.libPaths(c({local_sitelib}))
library(cn.mops)
library(rtracklayer)

count_drs <- getReadCountsFromBAM(bam_files, sampleNames=sample_names, mode="{pairmode}")
prep_counts <- cn.mops(count_drs)
cnv_out <- calcIntegerCopyNumbers(prep_counts)
cnvs <- cnvr(cnv_out)
export.bed(cnvs, {out_file})
"""
