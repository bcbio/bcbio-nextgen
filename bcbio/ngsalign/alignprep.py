"""Prepare read inputs (fastq, gzipped fastq and BAM) for parallel NGS alignment.
"""
import collections
import copy
import os
import subprocess

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils, qcsummary
from bcbio.provenance import do

def create_inputs(data):
    """Index input reads and prepare groups of reads to process concurrently.

    Allows parallelization of alignment beyond processors available on a single
    machine. Uses gbzip and grabix to prepare an indexed fastq file.
    """
    # skip skipping on samples without input files
    if data["files"][0] is None or data["algorithm"].get("align_split_size") is None:
        return [[data]]
    ready_files = _prep_grabix_indexes(data["files"], data["dirs"], data["config"])
    data["files"] = ready_files
    splits = _find_read_splits(ready_files[0], data["algorithm"]["align_split_size"])
    if len(splits) == 1:
        return [[data]]
    else:
        out = []
        for split in splits:
            cur_data = copy.deepcopy(data)
            cur_data["align_split"] = split
            out.append([cur_data])
        return out

def split_namedpipe_cl(in_file, data):
    """Create a commandline suitable for use as a named pipe with reads in a given region.
    """
    grabix = config_utils.get_program("grabix", data["config"])
    start, end = data["align_split"]
    return "<({grabix} grab {in_file} {start} {end})".format(**locals())

# ## configuration

def parallel_multiplier(items):
    """Determine if we will be parallelizing items during processing.
    """
    multiplier = 1
    for data in (x[0] for x in items):
        if data["algorithm"].get("align_split_size"):
            multiplier += 50
    return multiplier

# ## merge

def setup_combine(final_file, data):
    """Setup the data and outputs to allow merging data back together.
    """
    align_dir = os.path.dirname(final_file)
    base, ext = os.path.splitext(os.path.basename(final_file))
    start, end = data["align_split"]
    out_file = os.path.join(utils.safe_makedir(os.path.join(align_dir, "split")),
                            "%s-%s_%s%s" % (base, start, end, ext))
    data["combine"] = {"work_bam": {"out": final_file, "extras": []}}
    return out_file, data

def merge_split_alignments(samples, run_parallel):
    """Manage merging split alignments back into a final working BAM file.
    """
    ready = []
    file_key = "work_bam"
    to_merge = collections.defaultdict(list)
    for data in (xs[0] for xs in samples):
        if data.get("combine"):
            to_merge[data["combine"][file_key]["out"]].append(data)
        else:
            ready.append([data])
    ready_merge = []
    for mgroup in to_merge.itervalues():
        cur_data = mgroup[0]
        del cur_data["align_split"]
        for x in mgroup[1:]:
            cur_data["combine"][file_key]["extras"].append(x[file_key])
        ready_merge.append([cur_data])
    merged = run_parallel("delayed_bam_merge", ready_merge)
    return merged + ready

# ## determine file sections

def _find_read_splits(in_file, split_size):
    """Determine sections of fastq files to process in splits.

    Assumes a 4 line order to input files (name, read, name, quality).
    grabix is 1-based inclusive, so return coordinates in that format.
    """
    gbi_file = in_file + ".gbi"
    with open(gbi_file) as in_handle:
        in_handle.next() # throw away
        num_lines = int(in_handle.next().strip())
    assert num_lines % 4 == 0, "Expected lines to be multiple of 4"
    split_lines = split_size * 4
    chunks = []
    last = 1
    for chunki in range(num_lines // split_lines + min(1, num_lines % split_lines)):
        new = last + split_lines - 1
        chunks.append((last, min(new, num_lines - 1)))
        last = new
        if chunki > 0:
            last += 1
    return chunks

# ## bgzip and grabix

def _prep_grabix_indexes(in_files, dirs, config):
    if in_files[0].endswith(".bam") and in_files[1] is None:
        out = _bgzip_from_bam(in_files[0], dirs, config)
    else:
        out = [_bgzip_from_fastq(x, dirs, config) if x else None for x in in_files]
    [_grabix_index(x, config) for x in out if x]
    return out

def _bgzip_from_bam(bam_file, dirs, config):
    """Create bgzipped fastq files from an input BAM file.
    """
    # tools
    bedtools = config_utils.get_program("bedtools", config)
    bgzip = config_utils.get_program("bgzip", config)
    #bgzip = _get_bgzip_cmd(config)
    samtools = config_utils.get_program("samtools", config)
    resources = config_utils.get_resources("samtools", config)
    num_cores = config["algorithm"].get("num_cores", 1)
    max_mem = resources.get("memory", "1G")
    # files
    work_dir = utils.safe_makedir(os.path.join(dirs["work"], "align_prep"))
    out_file_1 = os.path.join(work_dir, "%s-1.fq.gz" % os.path.splitext(os.path.basename(bam_file))[0])
    if qcsummary.is_paired(bam_file):
        out_file_2 = out_file_1.replace("-1.fq.gz", "-2.fq.gz")
        fq2 = "-fq2 >(%s -c /dev/stdin > %s)" % (bgzip, out_file_2)
    else:
        out_file_2 = None
        fq2 = ""
    if not utils.file_exists(out_file_1):
        with file_transaction(out_file_1) as tx_out_file:
            fq1_bgzip_cmd = "%s -c /dev/stdin > %s" % (bgzip, tx_out_file)
            sortprefix = "%s-sort" % os.path.splitext(tx_out_file)[0]
            cmd = ("{samtools} sort -n -o -l 0 -@ {num_cores} -m {max_mem} {bam_file} {sortprefix} "
                   "| {bedtools} bamtofastq -i /dev/stdin "
                   "-fq >({fq1_bgzip_cmd}) {fq2}")
            do.run(cmd.format(**locals()), "BAM to bgzipped fastq")
    return [x for x in [out_file_1, out_file_2] if x is not None]

def _grabix_index(in_file, config):
    grabix = config_utils.get_program("grabix", config)
    gbi_file = in_file + ".gbi"
    if not utils.file_exists(gbi_file):
        do.run([grabix, "index", in_file], "Index input with grabix")
    return gbi_file

def _bgzip_from_fastq(in_file, dirs, config):
    """Prepare a bgzipped file from a fastq input, potentially gzipped (or bgzipped already).
    """
    grabix = config_utils.get_program("grabix", config)
    if in_file.endswith(".gz"):
        needs_bgzip, needs_gunzip = _check_gzipped_input(in_file, grabix)
    else:
        needs_bgzip, needs_gunzip = True, False
    if needs_bgzip or needs_gunzip:
        out_file = _bgzip_file(in_file, dirs, config, needs_bgzip, needs_gunzip)
    else:
        out_file = in_file
    return out_file

def _get_bgzip_cmd(config):
    """Retrieve command to use for bgzip, trying to use parallel pbgzip if available.
    """
    try:
        pbgzip = config_utils.get_program("pbgzip", config)
        return "%s -n %s " % (pbgzip, config["algorithm"].get("num_cores", 1))
    except config_utils.CmdNotFound:
        return config_utils.get_program("bgzip", config)

def _bgzip_file(in_file, dirs, config, needs_bgzip, needs_gunzip):
    """Handle bgzip of input file, potentially gunzipping an existing file.
    """
    work_dir = utils.safe_makedir(os.path.join(dirs["work"], "align_prep"))
    out_file = os.path.join(work_dir, os.path.basename(in_file) +
                            (".gz" if not in_file.endswith(".gz") else ""))
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            assert needs_bgzip
            bgzip = _get_bgzip_cmd(config)
            if needs_gunzip:
                gunzip_cmd = "gunzip -c {in_file} |".format(**locals())
                bgzip_in = "/dev/stdin"
            else:
                gunzip_cmd = ""
                bgzip_in = in_file
            do.run("{gunzip_cmd} {bgzip} -c {bgzip_in} > {tx_out_file}".format(**locals()),
                   "bgzip input file")
    return out_file

def _check_gzipped_input(in_file, grabix):
    """Determine if a gzipped input file is blocked gzip or standard.
    """
    is_bgzip = subprocess.check_output([grabix, "check", in_file])
    if is_bgzip.strip() == "yes":
        return False, False
    else:
        return True, True
