import os

from bcbio import bam, utils
from bcbio.provenance import do
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.log import logger
from bcbio.distributed.transaction import file_transaction


mapped_filter_query = (
    "not unmapped and "
    "not mate_is_unmapped and "
    "not secondary_alignment and "
    "not failed_quality_control")

def make_command(data, cmd, bam_file, bed_file=None,
                 depth_thresholds=None, max_cov=None, query=None, multicore=True):
    sambamba = config_utils.get_program("sambamba", data["config"], default="sambamba")
    num_cores = dd.get_cores(data) if multicore else 1
    target = (" -L " + bed_file) if bed_file else ""
    thresholds = "".join([" -T" + str(d) for d in (depth_thresholds or [])])
    maxcov = (" -C " + str(max_cov)) if max_cov else ""
    if query is None:
        query = mapped_filter_query + " and not duplicate"
    return ("{sambamba} {cmd} -t {num_cores} {bam_file} "
            "{target} {thresholds} {maxcov} -F \"{query}\"").format(**locals())

def index(data, bam_fpath):
    cmdl = make_command(data, "index", bam_fpath)
    indexed_bam = bam_fpath + ".bai"
    if not utils.file_uptodate(indexed_bam, bam_fpath):
        do.run(cmdl, "Indexing BAM file using sambamba")
        if not utils.file_exists(indexed_bam):
            logger.error("Cannot index BAM file " + bam_fpath + " using sambamba.")
            return None
    return indexed_bam

def work_dir(data):
    return utils.safe_makedir(os.path.join(dd.get_work_dir(data), "coverage", dd.get_sample_name(data), "sambamba"))

def _count_in_bam(data, bam_file, query, keep_dups=True, bed_file=None, target_name=None):
    if not keep_dups:
        if query:
            query += " and not duplicate"
        else:
            query = "not duplicate"
    cmd_id = "num_" + (query.replace(" ", "_") or "reads")
    if bed_file is not None:
        target_name = target_name or os.path.basename(bed_file)
        cmd_id += "_on_" + target_name
    output_file = os.path.join(work_dir(data), cmd_id)

    if not utils.file_uptodate(output_file, bam_file):
        index(data, bam_file)
        with file_transaction(data, output_file) as tx_out_file:
            cmdline = (make_command(data, "view -c", bam_file, bed_file, query=query, multicore=False)
                       + " > " + tx_out_file)
            do.run(cmdline, "Counting " + query + " for " + bam_file + ((" on " + target_name) if target_name else ""))

    with open(output_file) as f:
        return int(f.read().strip())

def number_of_reads(data, bam_file, keep_dups=True):
    return _count_in_bam(data, bam_file, "not secondary_alignment", keep_dups)

def number_of_mapped_reads(data, bam_file, keep_dups=True, bed_file=None, target_name=None):
    return _count_in_bam(data, bam_file, mapped_filter_query,
                         keep_dups=keep_dups, bed_file=bed_file, target_name=target_name)
