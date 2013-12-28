"""Clean an input BAM file to work with downstream pipelines.

GATK and Picard based pipelines have specific requirements for
chromosome order, run group information and other BAM formatting.
This provides a pipeline to prepare and resort an input.
"""
import os


from bcbio import bam, broad, utils
from bcbio.distributed.transaction import file_transaction

def picard_prep(in_bam, names, ref_file, dirs, config):
    """Prepare input BAM using Picard and GATK cleaning tools.

    - ReorderSam to reorder file to reference
    - AddOrReplaceReadGroups to add read group information and coordinate sort
    - PrintReads to filters to remove problem records:
      - filterMBQ to remove reads with mismatching bases and base qualities
    """
    runner = broad.runner_from_config(config)
    work_dir = utils.safe_makedir(os.path.join(dirs["work"], "bamclean", names["sample"]))
    reorder_bam = os.path.join(work_dir, "%s-reorder.bam" %
                               os.path.splitext(os.path.basename(in_bam))[0])
    reorder_bam = runner.run_fn("picard_reorder", in_bam, ref_file, reorder_bam)
    rg_bam = _fix_rgs_and_sort(reorder_bam, names, ref_file, runner)
    return _filter_bad_reads(rg_bam, ref_file, runner, config)

def _fix_rgs_and_sort(in_bam, names, ref_file, runner):
    """Fix input read groups if missing and coordinate sort file.
    """
    good_rg_bam = runner.run_fn("picard_fix_rgs", in_bam, names)
    return runner.run_fn("picard_sort", good_rg_bam, "coordinate")

def _filter_bad_reads(in_bam, ref_file, runner, config):
    """Use GATK filter to remove problem reads which choke GATK and Picard.
    """
    bam.index(in_bam, config)
    out_file = "%s-gatkfilter.bam" % os.path.splitext(in_bam)[0]
    if not utils.file_exists(out_file):
        with utils.curdir_tmpdir() as tmp_dir:
            with file_transaction(out_file) as tx_out_file:
                params = ["-T", "PrintReads",
                          "-R", ref_file,
                          "-I", in_bam,
                          "--out", tx_out_file,
                          "--filter_mismatching_base_and_quals"]
                runner.run_gatk(params, tmp_dir)
    return out_file
