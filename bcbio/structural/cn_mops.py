"""Copy number detection using read counts, with cn.mops.

http://www.bioconductor.org/packages/release/bioc/html/cn.mops.html
"""
from contextlib import closing
import os
import shutil
import subprocess

import pybedtools
import pysam

from bcbio import bam, install, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
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
    # XXX Can parallelize cn.mops with snow but currently causing
    # errors. Could also parallelize by chromosome.
    #num_cores = data["config"]["algorithm"].get("num_cores", 1)
    work_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural", names[0],
                                               "cn_mops"))
    with closing(pysam.Samfile(work_bams[0], "rb")) as pysam_work_bam:
        refs = pysam_work_bam.references
        # restrict to autosomes due to errors on X,Y and other chromosomes
        try:
            refs = refs[:refs.index("X")]
        except ValueError:
            refs = refs[:refs.index("chrX")]
        out_files = [_run_on_chrom(chrom, work_bams, names, work_dir, data)
                     for chrom in refs]
    out_file = _combine_out_files(out_files, work_bams[0], work_dir)
    out = []
    for data in items:
        if "sv" not in data:
            data["sv"] = {}
        data["sv"]["cnv"] = _prep_sample_cnvs(out_file, data)
        out.append(data)
    return out

def _combine_out_files(chr_files, base_bam, work_dir):
    """Concatenate all CNV calls into a single file.
    """
    out_file = os.path.join(work_dir, "%s-cnv.bed" % (os.path.splitext(os.path.basename(base_bam))[0]))
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                for chr_file in chr_files:
                    with open(chr_file) as in_handle:
                        is_empty = in_handle.readline().startswith("track name=empty")
                    if not is_empty:
                        with open(chr_file) as in_handle:
                            shutil.copyfileobj(in_handle, out_handle)
    return out_file

def _run_on_chrom(chrom, work_bams, names, work_dir, data):
    """Run cn.mops on work BAMs for a specific chromosome.
    """
    local_sitelib = os.path.join(install.get_defaults().get("tooldir", "/usr/local"),
                                 "lib", "R", "site-library")
    out_file = os.path.join(work_dir, "%s-%s-cnv.bed" % (os.path.splitext(os.path.basename(work_bams[0]))[0],
                                                         chrom))
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            rcode = "%s-run.R" % os.path.splitext(tx_out_file)[0]
            with open(rcode, "w") as out_handle:
                out_handle.write(_script.format(bam_file_str=",".join(work_bams), names_str=",".join(names),
                                                chrom=chrom, out_file=tx_out_file, num_cores=0,
                                                pairmode="paired" if bam.is_paired(work_bams[0]) else "unpaired",
                                                local_sitelib=local_sitelib))
            rscript = config_utils.get_program("Rscript", data["config"])
            try:
                do.run([rscript, rcode], "cn.mops CNV detection", data, log_error=False)
            except subprocess.CalledProcessError, msg:
                # cn.mops errors out if no CNVs found. Just write an empty file.
                if str(msg).find("No CNV regions in result object. Rerun cn.mops with different parameters") >= 0:
                    with open(tx_out_file, "w") as out_handle:
                        out_handle.write('track name=empty description="No CNVs found\n"')
                else:
                    logger.exception()
                    raise
    return out_file

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
.libPaths(c("{local_sitelib}"))
library(cn.mops)
library(rtracklayer)

bam_files <- strsplit("{bam_file_str}", ",")[[1]]
sample_names <- strsplit("{names_str}", ",")[[1]]
count_drs <- getReadCountsFromBAM(bam_files, sampleNames=sample_names, mode="{pairmode}",
                                  refSeqName={chrom}, parallel={num_cores})
prep_counts <- cn.mops(count_drs, parallel={num_cores})

cnv_out <- calcIntegerCopyNumbers(prep_counts)
calc_cnvs <- cnvs(cnv_out)
strcn_to_cn <- function(x) {{
  as.integer(substring(x, 3, 20))}}
calc_cnvs$score <- strcn_to_cn(calc_cnvs$CN)
calc_cnvs$name <- calc_cnvs$sampleName
export.bed(calc_cnvs, "{out_file}")
"""
