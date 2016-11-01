import os

from bcbio import utils
from bcbio.provenance import do
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.log import logger
from bcbio.distributed.transaction import file_transaction


def make_command(data, cmd, bam_file, bed_file=None,
                 depth_thresholds=None, max_cov=None, query=None):
    sambamba = config_utils.get_program("sambamba", data["config"], default="sambamba")
    num_cores = dd.get_cores(data)
    target = (" -L " + bed_file) if bed_file else ""
    thresholds = "".join([" -T" + str(d) for d in (depth_thresholds or [])])
    maxcov = (" -C " + str(max_cov)) if max_cov else ""
    if query is None:
        query = "not failed_quality_control and not duplicate and not unmapped"
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
    return utils.safe_makedir(os.path.join(dd.get_work_dir(data), "sambamba", dd.get_sample_name(data)))

def _count_in_bam(data, bam_file, query, keep_dups=True, bed_file=None, target_name=None):
    if not keep_dups:
        query += " and not duplicate"
    cmd_name = "num_" + (query.replace(" ", "_") or "reads")
    if bed_file is not None:
        target_name = target_name or os.path.basename(bed_file)
        cmd_name += "_on_" + target_name
    output_file = os.path.join(work_dir(data), cmd_name)

    if not utils.file_uptodate(output_file, bam_file):
        index(data, bam_file)
        with file_transaction(data, output_file) as tx_out_file:
            cmdline = make_command(data, "view -c", bam_file, bed_file, query=query) + " > " + tx_out_file
            do.run(cmdline, "Counting " + query + " for " + bam_file + ((" on " + target_name) if target_name else ""))

    with open(output_file) as f:
        return int(f.read().strip())

def number_of_reads(data, bam_file, keep_dups=True):
    return _count_in_bam(data, bam_file, '', keep_dups)

def number_of_mapped_reads(data, bam_file, keep_dups=True):
    return _count_in_bam(data, bam_file, 'not unmapped', keep_dups)

def number_of_properly_paired_reads(data, bam_file, keep_dups=True):
    return _count_in_bam(data, bam_file, 'proper_pair', keep_dups)

def number_of_dup_reads(data, bam_file):
    return _count_in_bam(data, bam_file, 'not unmapped and duplicate')

def number_mapped_reads_on_target(data, bed_file, bam_file, keep_dups=True, target_name=None):
    return _count_in_bam(data, bam_file, 'not unmapped', keep_dups, bed_file=bed_file, target_name=target_name)