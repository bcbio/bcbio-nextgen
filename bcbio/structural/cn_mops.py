"""Copy number detection using read counts, with cn.mops.

http://www.bioconductor.org/packages/release/bioc/html/cn.mops.html
"""
from contextlib import closing
import os
import re
import shutil
import subprocess

import pysam
import toolz as tz

from bcbio import bam, install, utils
from bcbio.distributed.multi import run_multicore, zeromq_aware_logging
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import config_utils, shared
from bcbio.provenance import do
from bcbio.structural import shared as sshared
from bcbio.variation import bedutils, vcfutils

def run(items, background=None):
    """Detect copy number variations from batched set of samples using cn.mops.
    """
    if not background: background = []
    names = [tz.get_in(["rgnames", "sample"], x) for x in items + background]
    work_bams = [x["align_bam"] for x in items + background]
    if len(items + background) < 2:
        raise ValueError("cn.mops only works on batches with multiple samples")
    data = items[0]
    work_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural", names[0],
                                               "cn_mops"))
    parallel = {"type": "local", "cores": data["config"]["algorithm"].get("num_cores", 1),
                "progs": ["delly"]}
    with closing(pysam.Samfile(work_bams[0], "rb")) as pysam_work_bam:
        chroms = [None] if _get_regional_bed_file(items[0]) else pysam_work_bam.references
        out_files = run_multicore(_run_on_chrom, [(chrom, work_bams, names, work_dir, items)
                                                  for chrom in chroms],
                                  data["config"], parallel)
    out_file = _combine_out_files(out_files, work_dir, data)
    out = []
    for data in items:
        if "sv" not in data:
            data["sv"] = []
        data["sv"].append({"variantcaller": "cn_mops",
                           "vrn_file": _prep_sample_cnvs(out_file, data)})
        out.append(data)
    return out

def _combine_out_files(chr_files, work_dir, data):
    """Concatenate all CNV calls into a single file.
    """
    out_file = "%s.bed" % sshared.outname_from_inputs(chr_files)
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                for chr_file in chr_files:
                    with open(chr_file) as in_handle:
                        is_empty = in_handle.readline().startswith("track name=empty")
                    if not is_empty:
                        with open(chr_file) as in_handle:
                            shutil.copyfileobj(in_handle, out_handle)
    return out_file

def _prep_sample_cnvs(cnv_file, data):
    """Convert a multiple sample CNV file into a single BED file for a sample.

    Handles matching and fixing names where R converts numerical IDs (1234) into
    strings by adding an X (X1234), and converts other characters into '.'s.
    http://stat.ethz.ch/R-manual/R-devel/library/base/html/make.names.html
    """
    import pybedtools
    sample_name = tz.get_in(["rgnames", "sample"], data)
    def make_names(name):
        return re.sub("[^\w.]", '.', name)
    def matches_sample_name(feat):
        return (feat.name == sample_name or feat.name == "X%s" % sample_name or
                feat.name == make_names(sample_name))
    def update_sample_name(feat):
        feat.name = sample_name
        return feat
    sample_file = os.path.join(os.path.dirname(cnv_file), "%s-cnv.bed" % sample_name)
    if not utils.file_exists(sample_file):
        with file_transaction(data, sample_file) as tx_out_file:
            with shared.bedtools_tmpdir(data):
                pybedtools.BedTool(cnv_file).filter(matches_sample_name).each(update_sample_name).saveas(tx_out_file)
    return sample_file

@utils.map_wrap
@zeromq_aware_logging
def _run_on_chrom(chrom, work_bams, names, work_dir, items):
    """Run cn.mops on work BAMs for a specific chromosome.
    """
    local_sitelib = os.path.join(install.get_defaults().get("tooldir", "/usr/local"),
                                 "lib", "R", "site-library")
    batch = sshared.get_cur_batch(items)
    ext = "-%s-cnv" % batch if batch else "-cnv"
    out_file = os.path.join(work_dir, "%s%s-%s.bed" % (os.path.splitext(os.path.basename(work_bams[0]))[0],
                                                       ext, chrom if chrom else "all"))
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            rcode = "%s-run.R" % os.path.splitext(out_file)[0]
            with open(rcode, "w") as out_handle:
                out_handle.write(_script.format(prep_str=_prep_load_script(work_bams, names, chrom, items),
                                                out_file=tx_out_file,
                                                local_sitelib=local_sitelib))
            rscript = utils.Rscript_cmd()
            try:
                do.run([rscript, rcode], "cn.mops CNV detection", items[0], log_error=False)
            except subprocess.CalledProcessError, msg:
                # cn.mops errors out if no CNVs found. Just write an empty file.
                if _allowed_cnmops_errorstates(str(msg)):
                    with open(tx_out_file, "w") as out_handle:
                        out_handle.write('track name=empty description="No CNVs found"\n')
                else:
                    logger.exception()
                    raise
    return [out_file]

def _allowed_cnmops_errorstates(msg):
    return (msg.find("No CNV regions in result object. Rerun cn.mops with different parameters") >= 0
            or msg.find("Normalization might not be applicable for this small number of segments") >= 0
            or msg.find("Error in if (is.finite(mv2m)) { : argument is of length zero") >= 0
            or msg.find("Some normalization factors are zero") >= 0)

def _prep_load_script(work_bams, names, chrom, items):
    if not chrom: chrom = ""
    pairmode = "paired" if bam.is_paired(work_bams[0]) else "unpaired"
    if len(items) == 2 and vcfutils.get_paired_phenotype(items[0]):
        load_script = _paired_load_script
    else:
        load_script = _population_load_script
    return load_script(work_bams, names, chrom, pairmode, items)

def _get_regional_bed_file(data):
    """If we are running a non-genome analysis, pull the regional file for analysis.
    """
    variant_regions = bedutils.merge_overlaps(tz.get_in(["config", "algorithm", "variant_regions"], data),
                                              data)
    is_genome = data["config"]["algorithm"].get("coverage_interval", "exome").lower() in ["genome"]
    if variant_regions and utils.file_exists(variant_regions) and not is_genome:
        return variant_regions

def _population_load_script(work_bams, names, chrom, pairmode, items):
    """Prepare BAMs for assessing CNVs in a population.
    """
    bed_file = _get_regional_bed_file(items[0])
    if bed_file:
        return _population_prep_targeted.format(bam_file_str=",".join(work_bams), names_str=",".join(names),
                                                chrom=chrom, num_cores=0, pairmode=pairmode, bed_file=bed_file)
    else:
        return _population_prep.format(bam_file_str=",".join(work_bams), names_str=",".join(names),
                                       chrom=chrom, num_cores=0, pairmode=pairmode)

def _paired_load_script(work_bams, names, chrom, pairmode, items):
    """Prepare BAMs for assessing CNVs in a paired tumor/normal setup.
    """
    paired = vcfutils.get_paired_bams(work_bams, items)
    bed_file = _get_regional_bed_file(items[0])
    if bed_file:
        return _paired_prep_targeted.format(case_file=paired.tumor_bam, case_name=paired.tumor_name,
                                            ctrl_file=paired.normal_bam, ctrl_name=paired.normal_name,
                                            num_cores=0, chrom=chrom, pairmode=pairmode, bed_file=bed_file)
    else:
        return _paired_prep.format(case_file=paired.tumor_bam, case_name=paired.tumor_name,
                                   ctrl_file=paired.normal_bam, ctrl_name=paired.normal_name,
                                   num_cores=0, chrom=chrom, pairmode=pairmode)

_script = """
.libPaths(c("{local_sitelib}"))
library(cn.mops)
library(rtracklayer)

{prep_str}

calc_cnvs <- cnvs(cnv_out)
strcn_to_cn <- function(x) {{
  as.numeric(substring(x, 3, 20))}}
calc_cnvs$score <- strcn_to_cn(calc_cnvs$CN)
calc_cnvs$name <- calc_cnvs$sampleName
export.bed(calc_cnvs, "{out_file}")
"""

_population_prep = """
bam_files <- strsplit("{bam_file_str}", ",")[[1]]
sample_names <- strsplit("{names_str}", ",")[[1]]
count_drs <- getReadCountsFromBAM(bam_files, sampleNames=sample_names, mode="{pairmode}",
                                  refSeqName="{chrom}", parallel={num_cores})
prep_counts <- cn.mops(count_drs, parallel={num_cores})
cnv_out <- calcIntegerCopyNumbers(prep_counts)
"""

_paired_prep = """
case_count <- getReadCountsFromBAM(c("{case_file}"), sampleNames=c("{case_name}"), mode="{pairmode}",
                                   refSeqName="{chrom}", parallel={num_cores})
ctrl_count <- getReadCountsFromBAM(c("{ctrl_file}"), sampleNames=c("{ctrl_name}"), mode="{pairmode}",
                                   refSeqName="{chrom}", parallel={num_cores},
                                   WL=width(case_count)[[1]])
prep_counts <- referencecn.mops(case_count, ctrl_count, parallel={num_cores})
cnv_out <- calcIntegerCopyNumbers(prep_counts)
"""

_population_prep_targeted = """
bam_files <- strsplit("{bam_file_str}", ",")[[1]]
sample_names <- strsplit("{names_str}", ",")[[1]]
my_gr <- import.bed(c("{bed_file}"), trackLine=FALSE, asRangedData=FALSE)
if ("{chrom}" != "") my_gr = subset(my_gr, seqnames(my_gr) == "{chrom}")
if (length(my_gr) < 1) stop("No CNV regions in result object. Rerun cn.mops with different parameters!")
count_drs <- getSegmentReadCountsFromBAM(bam_files, sampleNames=sample_names, mode="{pairmode}",
                                         GR=my_gr, parallel={num_cores})
prep_counts <- cn.mops(count_drs, parallel={num_cores})
cnv_out <- calcIntegerCopyNumbers(prep_counts)
"""

_paired_prep_targeted = """
my_gr <- import.bed(c("{bed_file}"), trackLine=FALSE, asRangedData=FALSE)
if ("{chrom}" != "") my_gr = subset(my_gr, seqnames(my_gr) == "{chrom}")
if (length(my_gr) < 1) stop("No CNV regions in result object. Rerun cn.mops with different parameters!")
case_count <- getSegmentReadCountsFromBAM(c("{case_file}"), GR=my_gr,
                                          sampleNames=c("{case_name}"),
                                          mode="{pairmode}", parallel={num_cores})
ctrl_count <- getSegmentReadCountsFromBAM(c("{ctrl_file}"), GR=my_gr,
                                          sampleNames=c("{case_name}"),
                                          mode="{pairmode}", parallel={num_cores})
prep_counts <- referencecn.mops(case_count, ctrl_count, parallel={num_cores})
cnv_out <- calcIntegerCopyNumbers(prep_counts)
"""
