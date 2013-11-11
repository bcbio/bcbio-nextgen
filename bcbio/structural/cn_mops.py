"""Copy number detection using read counts, with cn.mops.

http://www.bioconductor.org/packages/release/bioc/html/cn.mops.html
"""
import os

import pybedtools

from bcbio import bam, install, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.provenance import do

def run(items):
    """Detect copy number variations from batched set of samples using cn.mops.
    """
    names = [x["name"][-1] for x in items]
    work_bams = [x["work_bam"] for x in items]
    if len(items) == 1:
        raise ValueError("cn.mops only works on batches with multiple samples")
    data = items[0]
    num_cores = data["config"]["algorithm"].get("num_cores", 1)
    work_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural", names[0],
                                               "cn_mops"))
    local_sitelib = os.path.join(install.get_defaults().get("tooldir", "/usr/local"),
                                 "lib", "R", "site-library")
    out_file = os.path.join(work_dir, "%s-cnv.bed" % (os.path.splitext(os.path.basename(work_bams[0]))[0]))
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            rcode = "%s-run.R" % os.path.splitext(tx_out_file)[0]
            with open(rcode, "w") as out_handle:
                out_handle.write(_script.format(bam_file_str=",".join(work_bams), names_str=",".join(names),
                                                out_file=tx_out_file, num_cores=0 if num_cores == 1 else num_cores,
                                                pairmode="paired" if bam.is_paired(work_bams[0]) else "unpaired",
                                                local_sitelib=local_sitelib))
            rscript = config_utils.get_program("Rscript", data["config"])
            do.run([rscript, rcode], "cn.mops CNV detection", data)
    out = []
    for data in items:
        if "sv" not in data:
            data["sv"] = {}
        data["sv"]["cnv"] = _prep_sample_cnvs(out_file, data)
        out.append(data)
    return out

def _prep_sample_cnvs(cnv_file, data):
    """Convert a multiple sample CNV file into
    """
    sample_name = data["name"][-1]
    sample_file = os.path.join(os.path.dirname(cnv_file), "%s-cnv.bed" % sample_name)
    if not utils.file_exists(sample_file):
        with file_transaction(sample_file) as tx_out_file:
            pybedtools.BedTool(cnv_file).filter(lambda x: x.name == sample_name).saveas(tx_out_file)
    return sample_file

_script = """
bam_files <- strsplit("{bam_file_str}", ",")[[1]]
sample_names <- strsplit("{names_str}", ",")[[1]]

.libPaths(c("{local_sitelib}"))
library(cn.mops)
library(rtracklayer)

count_drs <- getReadCountsFromBAM(bam_files, sampleNames=sample_names, mode="{pairmode}",
                                  parallel={num_cores})
prep_counts <- cn.mops(count_drs, parallel={num_cores})
cnv_out <- calcIntegerCopyNumbers(prep_counts)
calc_cnvs <- cnvs(cnv_out)
strcn_to_cn <- function(x) {
  as.integer(substring(x, 3, 20))}
calc_cnvs$score <- strcn_to_cn(calc_cnvs$CN)
calc_cnvs$name <- calc_cnvs$sampleName
export.bed(calc_cnvs, "{out_file}")
"""
