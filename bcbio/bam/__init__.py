"""Functionality to query and extract information from aligned BAM files.
"""
import contextlib
import os
import subprocess

import pysam

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.provenance import do

def is_paired(bam_file):
    """Determine if a BAM file has paired reads.
    """
    with contextlib.closing(pysam.Samfile(bam_file, "rb")) as in_pysam:
        for read in in_pysam:
            return read.is_paired

def index(in_bam, config):
    """Index a BAM file, skipping if index present.

    Centralizes BAM indexing providing ability to switch indexing approaches.
    """
    assert is_bam(in_bam), "%s in not a BAM file" % in_bam
    index_file = "%s.bai" % in_bam
    alt_index_file = "%s.bai" % os.path.splitext(in_bam)[0]
    if (not utils.file_uptodate(index_file, in_bam) and
          not utils.file_uptodate(alt_index_file, in_bam)):
        sambamba = _get_sambamba(config)
        samtools = config_utils.get_program("samtools", config)
        num_cores = config["algorithm"].get("num_cores", 1)
        with file_transaction(index_file) as tx_index_file:
            samtools_cmd = "{samtools} index {in_bam} {tx_index_file}"
            if sambamba:
                cmd = "{sambamba} index -t {num_cores} {in_bam} {tx_index_file}"
            else:
                cmd = samtools_cmd
            # sambamba has intermittent multicore failures. Allow
            # retries with single core
            try:
                do.run(cmd.format(**locals()), "Index BAM file: %s" % os.path.basename(in_bam),
                       log_error=False)
            except:
                do.run(samtools_cmd.format(**locals()),
                       "Index BAM file (single core): %s" % os.path.basename(in_bam))
    return index_file if utils.file_uptodate(index_file, in_bam) else alt_index_file

def get_downsample_pct(runner, in_bam, target_counts):
    """Retrieve percentage of file to downsample to get to target counts.
    """
    total = sum(x.aligned for x in runner.run_fn("picard_idxstats", in_bam))
    with contextlib.closing(pysam.Samfile(in_bam, "rb")) as work_bam:
        n_rgs = max(1, len(work_bam.header["RG"]))
    rg_target = n_rgs * target_counts
    if total > rg_target:
        return float(rg_target) / float(total)

def downsample(in_bam, data, target_counts):
    """Downsample a BAM file to the specified number of target counts.
    """
    broad_runner = broad.runner_from_config(data["config"])
    index(in_bam, data["config"])
    ds_pct = get_downsample_pct(broad_runner, in_bam, target_counts)
    if ds_pct:
        out_file = "%s-downsample%s" % os.path.splitext(in_bam)
        if not utils.file_exists(out_file):
            with file_transaction(out_file) as tx_out_file:
                sambamba = config_utils.get_program("sambamba", data["config"])
                num_cores = data["config"]["algorithm"].get("num_cores", 1)
                cmd = ("{sambamba} view -t {num_cores} -f bam -o {tx_out_file} "
                       "--subsample={ds_pct:.3} --subsampling-seed=42 {in_bam}")
                do.run(cmd.format(**locals()), "Downsample BAM file: %s" % os.path.basename(in_bam))
        return out_file

def open_samfile(in_file):
    if is_bam(in_file):
        return pysam.Samfile(in_file, "rb")
    elif is_sam(in_file):
        return pysam.Samfile(in_file, "r")
    else:
        raise IOError("in_file must be either a BAM file or SAM file. Is the "
                      "extension .sam or .bam?")

def is_bam(in_file):
    _, ext = os.path.splitext(in_file)
    if ext == ".bam":
        return True
    else:
        return False


def is_sam(in_file):
    _, ext = os.path.splitext(in_file)
    if ext == ".sam":
        return True
    else:
        return False

def count(in_bam, config=None):
    """
    return the counts in a BAM file
    """
    if not config:
        config = {}
    sambamba = _get_sambamba(config)
    if sambamba:
        cmd = ("{sambamba} view -c {in_bam}").format(**locals())
    else:
        samtools = config_utils.get_program("samtools", config)
        cmd = ("{samtools} view -c {in_bam}").format(**locals())
    out = subprocess.check_output(cmd, shell=True)
    return int(out)

def sam_to_bam(in_sam, config):
    if is_bam(in_sam):
        return in_sam

    assert is_sam(in_sam), "%s is not a SAM file" % in_sam
    out_file = os.path.splitext(in_sam)[0] + ".bam"
    if utils.file_exists(out_file):
        return out_file

    samtools = config_utils.get_program("samtools", config)
    num_cores = config["algorithm"].get("num_cores", 1)
    with file_transaction(out_file) as tx_out_file:
        cmd = "{samtools} view -@ {num_cores} -h -S -b {in_sam} -o {tx_out_file}"
        do.run(cmd.format(**locals()),
               ("Convert SAM to BAM (%s cores): %s to %s"
                % (str(num_cores), in_sam, out_file)))
    return out_file

def bam_to_sam(in_file, config):
    if is_sam(in_file):
        return in_file

    assert is_bam(in_file), "%s is not a BAM file" % in_file
    out_file = os.path.splitext(in_file)[0] + ".sam"
    if utils.file_exists(out_file):
        return out_file

    samtools = config_utils.get_program("samtools", config)
    num_cores = config["algorithm"].get("num_cores", 1)
    with file_transaction(out_file) as tx_out_file:
        cmd = "{samtools} view -@ {num_cores} -h {in_file} -o {tx_out_file}"
        do.run(cmd.format(**locals()),
               ("Convert BAM to SAM (%s cores): %s to %s"
                % (str(num_cores), in_file, out_file)))
    return out_file


def merge(bamfiles, out_bam, config):
    assert all(map(is_bam, bamfiles)), ("Not all of the files to merge are not BAM "
                                        "files: %s " % (bamfiles))
    assert all(map(utils.file_exists, bamfiles)), ("Not all of the files to merge "
                                             "exist: %s" % (bamfiles))
    sambamba = _get_sambamba(config)
    sambamba = None
    samtools = config_utils.get_program("samtools", config)
    num_cores = config["algorithm"].get("num_cores", 1)
    with file_transaction(out_bam) as tx_out_bam:
        if sambamba:
            cmd = "{sambamba} merge -t {num_cores} {tx_out_bam} " + " ".join(bamfiles)
        else:
            cmd = "{samtools} merge -@ {num_cores} {tx_out_bam} " + " ".join(bamfiles)
        do.run(cmd.format(**locals()), "Merge %s into %s." % (bamfiles, out_bam))
    return out_bam


def sort(in_bam, config, order="coordinate"):
    """Sort a BAM file, skipping if already present.
    """
    assert is_bam(in_bam), "%s in not a BAM file" % in_bam
    if bam_already_sorted(in_bam, config, order):
        return in_bam

    sort_stem = _get_sort_stem(in_bam, order)
    sort_file = sort_stem + ".bam"
    if not utils.file_exists(sort_file):
        sambamba = _get_sambamba(config)
        sambamba = None
        samtools = config_utils.get_program("samtools", config)
        num_cores = config["algorithm"].get("num_cores", 1)
        with file_transaction(sort_file) as tx_sort_file:
            tx_sort_stem = os.path.splitext(tx_sort_file)[0]
            tx_dir = utils.safe_makedir(os.path.dirname(tx_sort_file))
            order_flag = "-n" if order is "queryname" else ""
            samtools_cmd = ("{samtools} sort {order_flag} "
                            "{in_bam} {tx_sort_stem}")
            if sambamba:
                cmd = ("{sambamba} sort -t {num_cores} {order_flag} "
                       "-o {tx_sort_file} --tmpdir={tx_dir} {in_bam}")
            else:
                cmd = samtools_cmd
            # sambamba has intermittent multicore failures. Allow
            # retries with single core
            try:
                do.run(cmd.format(**locals()),
                       "Sort BAM file (multi core, %s): %s to %s" %
                       (order, os.path.basename(in_bam),
                        os.path.basename(sort_file)), log_error=False)
            except:
                do.run(samtools_cmd.format(**locals()),
                       "Sort BAM file (single core, %s): %s" %
                       (order, os.path.basename(in_bam),
                        os.path.basename(sort_file)))
    return sort_file


def _get_sambamba(config):
    try:
        sambamba = config_utils.get_program("sambamba", config)
    except config_utils.CmdNotFound:
        sambamba = None
    return sambamba


def bam_already_sorted(in_bam, config, order):
    return order == _get_sort_order(in_bam, config)


def _get_sort_order(in_bam, config):
    with open_samfile(in_bam) as bam_handle:
        header = bam_handle.header
    return utils.get_in(header, ("HD", "SO"), None)

def _get_sort_stem(in_bam, order):
    SUFFIXES = {"coordinate": ".sorted", "queryname": ".nsorted"}
    sort_base = os.path.splitext(in_bam)[0]
    for suffix in SUFFIXES:
        sort_base = sort_base.split(suffix)[0]
    return sort_base + SUFFIXES[order]

def sample_name(in_bam):
    """
    get sample name from a BAM file as a work around for pysam 0.6 header
    parsing issues. This can get replaced by the pysam version in cortex.py
    when we move to Docker.
    """
    cmd = "samtools view -H {in_bam}".format(**locals())
    out, _ = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()
    name = None
    for line in out.split("\n"):
        if line.startswith("@RG"):
            name = [x.split(":")[1] for x in line.split() if x.split(":")[0] == "SM"]
    return name[0] if name else None
