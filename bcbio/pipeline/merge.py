"""Handle multiple samples present on a single flowcell

Merges samples located in multiple lanes on a flowcell. Unique sample names identify
items to combine within a group.
"""
import os
import shutil
import subprocess

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.pipeline import config_utils
from bcbio.provenance import do, system

def combine_fastq_files(in_files, work_dir, config):
    if len(in_files) == 1:
        return in_files[0]
    else:
        cur1, cur2 = in_files[0]
        out1 = os.path.join(work_dir, os.path.basename(cur1))
        out2 = os.path.join(work_dir, os.path.basename(cur2)) if cur2 else None
        if not os.path.exists(out1):
            with open(out1, "a") as out_handle:
                for (cur1, _) in in_files:
                    with open(cur1) as in_handle:
                        shutil.copyfileobj(in_handle, out_handle)
        if out2 and not os.path.exists(out2):
            with open(out2, "a") as out_handle:
                for (_, cur2) in in_files:
                    with open(cur2) as in_handle:
                        shutil.copyfileobj(in_handle, out_handle)
        for f1, f2 in in_files:
            utils.save_diskspace(f1, "fastq merged to %s" % out1, config)
            if f2:
                utils.save_diskspace(f2, "fastq merged to %s" % out2, config)
        return out1, out2

def merge_bam_files(bam_files, work_dir, config, out_file=None, batch=None):
    """Merge multiple BAM files from a sample into a single BAM for processing.

    Checks system open file limit and merges in batches if necessary to avoid
    file handle limits.
    """
    if len(bam_files) == 1:
        bam.index(bam_files[0], config)
        return bam_files[0]
    else:
        if out_file is None:
            out_file = os.path.join(work_dir, os.path.basename(sorted(bam_files)[0]))
        if batch is not None:
            base, ext = os.path.splitext(out_file)
            out_file = "%s-b%s%s" % (base, batch, ext)
        if not utils.file_exists(out_file):
            sambamba = config_utils.get_program("sambamba", config)
            samtools = config_utils.get_program("samtools", config)
            resources = config_utils.get_resources("samtools", config)
            num_cores = config["algorithm"].get("num_cores", 1)
            max_mem = config_utils.adjust_memory(resources.get("memory", "1G"),
                                                 2, "decrease").upper()
            # sambamba opens 4 handles per file, so try to guess a reasonable batch size
            batch_size = (system.open_file_limit() // 4) - 100
            if len(bam_files) > batch_size:
                bam_files = [merge_bam_files(xs, work_dir, config, out_file, i)
                             for i, xs in enumerate(utils.partition_all(batch_size, bam_files))]
            with tx_tmpdir(config) as tmpdir:
                with utils.chdir(tmpdir):
                    with file_transaction(config, out_file) as tx_out_file:
                        with file_transaction(config, "%s.list" % os.path.splitext(out_file)[0]) as tx_bam_file_list:
                            with open(tx_bam_file_list, "w") as out_handle:
                                for f in sorted(bam_files):
                                    out_handle.write("%s\n" % f)
                            if bam.bam_already_sorted(bam_files[0], config, "coordinate"):
                                cmd = _sambamba_merge(bam_files)
                            else:
                                assert config.get("mark_duplicates", True)
                                cmd = _biobambam_merge_dedup()
                            do.run(cmd.format(**locals()), "Merge bam files to %s" % os.path.basename(out_file),
                                   None)
            # Ensure timestamps are up to date on output file and index
            # Works around issues on systems with inconsistent times
            for ext in ["", ".bai"]:
                if os.path.exists(out_file + ext):
                    subprocess.check_call(["touch", out_file + ext])
            for b in bam_files:
                utils.save_diskspace(b, "BAM merged to %s" % out_file, config)
        bam.index(out_file, config)
        return out_file

def _biobambam_merge_dedup():
    """Combine query sorted BAM files, de-duplicate and sort. Handles split prepped files.
    """
    return ("bammerge IL={tx_bam_file_list} SO=queryname level=0 tmpfile={tx_out_file}-bammerge | "
            "bamsormadup threads={num_cores} tmpfile={tx_out_file}-bamsormaduptmp > {tx_out_file}")

def _sambamba_merge(bam_files):
    """Merge multiple BAM files with sambamba.
    """
    if len(bam_files) > system.open_file_limit():
        raise IOError("More files to merge (%s) than available open file descriptors (%s)\n"
                      "See documentation on tips for changing file limits:\n"
                      "https://bcbio-nextgen.readthedocs.org/en/latest/contents/"
                      "parallel.html#tuning-systems-for-scale"
                      % (len(bam_files), system.open_file_limit()))
    return "{sambamba} merge {tx_out_file} -t {num_cores} `cat {tx_bam_file_list}`"
