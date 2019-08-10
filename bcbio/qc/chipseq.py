"""Qc for chipseq pipeline"""
import os
import shutil
import glob

from bcbio.log import logger
from bcbio import utils
from bcbio.bam.readstats import number_of_mapped_reads
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.distributed.transaction import tx_tmpdir


supported = ["hg19", "hg38", "mm10", "mm9", "rn4", "ce6", "dm3"]

def run(bam_file, sample, out_dir):
    """Standard QC metrics for chipseq"""
    out = {}
    # if "rchipqc" in dd.get_tools_on(sample):
    #    out = chipqc(bam_file, sample, out_dir)

    peaks = sample.get("peaks_files", {}).get("main")
    if peaks:
        out.update(_reads_in_peaks(bam_file, peaks, sample))
    return out

def _reads_in_peaks(bam_file, peaks_file, sample):
    """Calculate number of reads in peaks"""
    if not peaks_file:
        return {}
    rip = number_of_mapped_reads(sample, bam_file, bed_file = peaks_file)
    return {"metrics": {"RiP": rip}}

def chipqc(bam_file, sample, out_dir):
    """Attempt code to run ChIPQC bioconductor packate in one sample"""
    sample_name = dd.get_sample_name(sample)
    logger.warning("ChIPQC is unstable right now, if it breaks, turn off the tool.")
    if utils.file_exists(out_dir):
        return _get_output(out_dir)
    with tx_tmpdir() as tmp_dir:
        rcode = _sample_template(sample, tmp_dir)
        if rcode:
            # local_sitelib = utils.R_sitelib()
            rscript = utils.Rscript_cmd()
            do.run([rscript, "--vanilla", rcode], "ChIPQC in %s" % sample_name, log_error=False)
            shutil.move(tmp_dir, out_dir)
    return _get_output(out_dir)

def _get_output(out_dir):
    files = glob.glob(out_dir)
    if files:
        return {'secondary': files}
    else:
        return {}

def _sample_template(sample, out_dir):
    """R code to get QC for one sample"""
    bam_fn = dd.get_work_bam(sample)
    genome = dd.get_genome_build(sample)
    if genome in supported:
        peaks = sample.get("peaks_files", []).get("main")
        if peaks:
            r_code = ("library(ChIPQC);\n"
                      "sample = ChIPQCsample(\"{bam_fn}\","
                      "\"{peaks}\", "
                      "annotation = \"{genome}\","
                      ");\n"
                      "ChIPQCreport(sample);\n")
            r_code_fn = os.path.join(out_dir, "chipqc.r")
            with open(r_code_fn, 'w') as inh:
                inh.write(r_code.format(**locals()))
            return r_code_fn

def _template():
    r_code = ("library(ChIPQC);\n"
              "metadata = read.csv(\"metadata.csv\");\n"
              "qc = ChIPQC(metadata, annotation = {genome},"
              "     chromosomes = {chr}, fragmentLength = {fl});\n"
              "ChIPQCreport(qc, facetBy=c(\"{facet1}\", \"{facet2}\"));\n")
