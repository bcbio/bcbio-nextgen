import os
import subprocess

from bcbio import bam, utils
from bcbio.provenance import do
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.distributed.transaction import file_transaction

pybedtools = utils.LazyImport("pybedtools")

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

def work_dir(data):
    return utils.safe_makedir(os.path.join(dd.get_work_dir(data), "coverage", dd.get_sample_name(data), "sambamba"))

def number_of_mapped_reads(data, bam_file, keep_dups=True, bed_file=None, target_name=None):
    """Count mapped reads, allow adjustment for duplicates and BED regions.

    Since samtools view does not use indexes for BED files
    (https://github.com/samtools/samtools/issues/88)
    we loop over regions in a BED file and add the counts together.
    """
    # Flag explainer https://broadinstitute.github.io/picard/explain-flags.html
    if keep_dups:
        query_name = mapped_filter_query
        flag = 780  # not (read unmapped or mate unmapped or fails QC or secondary alignment)
    else:
        query_name = mapped_filter_query + " and not duplicate"
        flag = 1804  # as above plus not duplicate
    cmd_id = "num_" + query_name.replace(" ", "_")
    if bed_file is not None:
        target_name = target_name or os.path.basename(bed_file)
        cmd_id += "_on_" + target_name
    output_file = os.path.join(work_dir(data), cmd_id)
    if not utils.file_uptodate(output_file, bam_file):
        bam.index(bam_file, data["config"], check_timestamp=False)
        num_cores = dd.get_num_cores(data)
        with file_transaction(data, output_file) as tx_out_file:
            count = 0
            cmd = "samtools view -c -F {flag} -@ {num_cores} {bam_file}{region}"
            for r in (pybedtools.BedTool(bed_file) if bed_file else [None]):
                # Covert to samtools region (these are 1-based, BED is 0-based)
                region = " %s:%s-%s" % (r.chrom, r.start + 1, r.end) if r else ""
                count += int(subprocess.check_output(cmd.format(**locals()), shell=True))
            with open(tx_out_file, "w") as out_handle:
                out_handle.write(str(count))
    with open(output_file) as f:
        return int(f.read().strip())
