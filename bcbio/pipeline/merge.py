"""Handle multiple samples present on a single flowcell

Merges samples located in multiple lanes on a flowcell. Unique sample names identify
items to combine within a group.
"""
import os
import shutil

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
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
        return bam_files[0]
    else:
        if out_file is None:
            out_file = os.path.join(work_dir, os.path.basename(sorted(bam_files)[0]))
        if batch is not None:
            base, ext = os.path.splitext(out_file)
            out_file = "%s-b%s%s" % (base, batch, ext)
        if not utils.file_exists(out_file) or not utils.file_exists(out_file + ".bai"):
            bamtools = config_utils.get_program("bamtools", config)
            samtools = config_utils.get_program("samtools", config)
            resources = config_utils.get_resources("samtools", config)
            num_cores = config["algorithm"].get("num_cores", 1)
            max_mem = resources.get("memory", "1G")
            batch_size = system.open_file_limit() - 100
            if len(bam_files) > batch_size:
                bam_files = [merge_bam_files(xs, work_dir, config, out_file, i)
                             for i, xs in enumerate(utils.partition_all(batch_size, bam_files))]
            with utils.curdir_tmpdir() as tmpdir:
                with utils.chdir(tmpdir):
                    merge_cl = _bamtools_merge(bam_files)
                    with file_transaction(out_file) as tx_out_file:
                        tx_out_prefix = os.path.splitext(tx_out_file)[0]
                        with utils.tmpfile(dir=work_dir, prefix="bammergelist") as bam_file_list:
                            bam_file_list = "%s.list" % os.path.splitext(out_file)[0]
                            with open(bam_file_list, "w") as out_handle:
                                for f in sorted(bam_files):
                                    out_handle.write("%s\n" % f)
                            cmd = (merge_cl + " | "
                                   "{samtools} sort -@ {num_cores} -m {max_mem} - {tx_out_prefix}")
                            do.run(cmd.format(**locals()), "Merge bam files", None)
            for b in bam_files:
                utils.save_diskspace(b, "BAM merged to %s" % out_file, config)
        bam.index(out_file, config)
        return out_file

def _samtools_cat(bam_files, tmpdir):
    """Concatenate multiple BAM files together with samtools.
    Creates short paths to shorten the commandline.
    """
    short_bams = []
    for i, bam_file in enumerate(bam_files):
        short_bam = os.path.join(tmpdir, "%s.bam" % i)
        os.symlink(bam_file, short_bam)
        short_bams.append(short_bam)
    return "{samtools} cat " + " ".join(os.path.relpath(b) for b in short_bams)

def _bamtools_merge(bam_files):
    """Use bamtools to merge multiple BAM files, requires a list from disk.
    """
    if len(bam_files) > system.open_file_limit():
        raise IOError("More files to merge (%s) than available open file descriptors (%s)\n"
                      "See documentation on tips for changing file limits:\n"
                      "https://bcbio-nextgen.readthedocs.org/en/latest/contents/"
                      "parallel.html#tuning-systems-for-scale"
                      % (len(bam_files), system.open_file_limit()))
    return "{bamtools} merge -list {bam_file_list}"
