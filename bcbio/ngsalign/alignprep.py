"""Prepare read inputs (fastq, gzipped fastq and BAM) for parallel NGS alignment.
"""
import collections
import copy
import os
import subprocess

import toolz as tz
try:
    import pybedtools
except ImportError:
    pybedtools = None

from bcbio import bam, utils
from bcbio.bam import cram
from bcbio.log import logger
from bcbio.distributed.multi import run_multicore, zeromq_aware_logging
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils, tools
from bcbio.provenance import do

def create_inputs(data):
    """Index input reads and prepare groups of reads to process concurrently.

    Allows parallelization of alignment beyond processors available on a single
    machine. Uses gbzip and grabix to prepare an indexed fastq file.
    """
    # CRAM files must be converted to bgzipped fastq
    if not ("files" in data and _is_cram_input(data["files"])):
        # skip indexing on samples without input files or not doing alignment
        if ("files" not in data or data["files"][0] is None or
              data["config"]["algorithm"].get("align_split_size") is None
              or not data["config"]["algorithm"].get("aligner")):
            return [[data]]
    ready_files = _prep_grabix_indexes(data["files"], data["dirs"], data)
    data["files"] = ready_files
    # bgzip preparation takes care of converting illumina into sanger format
    data["config"]["algorithm"]["quality_format"] = "standard"
    if tz.get_in(["config", "algorithm", "align_split_size"], data):
        splits = _find_read_splits(ready_files[0], data["config"]["algorithm"]["align_split_size"])
    else:
        splits = [None]
    if len(splits) == 1:
        return [[data]]
    else:
        out = []
        for split in splits:
            cur_data = copy.deepcopy(data)
            cur_data["align_split"] = list(split)
            out.append([cur_data])
        return out

def split_namedpipe_cl(in_file, data):
    """Create a commandline suitable for use as a named pipe with reads in a given region.
    """
    grabix = config_utils.get_program("grabix", data["config"])
    start, end = data["align_split"]
    return "<({grabix} grab {in_file} {start} {end})".format(**locals())

def fastq_convert_pipe_cl(in_file, data):
    """Create an anonymous pipe converting Illumina 1.3-1.7 to Sanger.

    Uses seqtk: https://github.com/lh3/seqt
    """
    seqtk = config_utils.get_program("seqtk", data["config"])
    return "<({seqtk} seq -Q64 -V {in_file})".format(**locals())

# ## configuration

def parallel_multiplier(items):
    """Determine if we will be parallelizing items during processing.
    """
    multiplier = 1
    for data in (x[0] for x in items):
        if data["config"]["algorithm"].get("align_split_size"):
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
        in_handle.next()  # throw away
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

def _is_bam_input(in_files):
    return in_files[0].endswith(".bam") and (len(in_files) == 1 or in_files[1] is None)

def _is_cram_input(in_files):
    return in_files[0].endswith(".cram") and (len(in_files) == 1 or in_files[1] is None)

def _prep_grabix_indexes(in_files, dirs, data):
    if _is_bam_input(in_files):
        out = _bgzip_from_bam(in_files[0], dirs, data["config"])
    elif _is_cram_input(in_files):
        out = _bgzip_from_cram(in_files[0], dirs, data)
    else:
        out = run_multicore(_bgzip_from_fastq,
                            [[{"in_file": x, "dirs": dirs, "config": data["config"]}] for x in in_files if x],
                            data["config"])
    items = [[{"bgzip_file": x, "config": copy.deepcopy(data["config"])}] for x in out if x]
    run_multicore(_grabix_index, items, data["config"])
    return out

def _bgzip_from_cram(cram_file, dirs, data):
    """Create bgzipped fastq files from an input CRAM file in regions of interest.

    Returns a list with a single file, for single end CRAM files, or two
    files for paired end input.
    """
    region_file = (tz.get_in(["config", "algorithm", "variant_regions"], data)
                   if tz.get_in(["config", "algorithm", "coverage_interval"], data) in ["regional", "exome"]
                   else None)
    if region_file:
        regions = ["%s:%s-%s" % tuple(r) for r in pybedtools.BedTool(region_file)]
    else:
        regions = [None]
    work_dir = utils.safe_makedir(os.path.join(dirs["work"], "align_prep"))
    out_s, out_p1, out_p2 = [os.path.join(work_dir, "%s-%s.fq.gz" %
                                          (utils.splitext_plus(os.path.basename(cram_file))[0], fext))
                             for fext in ["s1", "p1", "p2"]]
    if not utils.file_exists(out_s) and not utils.file_exists(out_p1):
        cram.index(cram_file)
        fastqs = _cram_to_fastq_regions(regions, cram_file, dirs, data)
        if len(fastqs[0]) == 1:
            with file_transaction(out_s) as tx_out_file:
                _merge_and_bgzip([xs[0] for xs in fastqs], tx_out_file, out_s)
        else:
            for i, out_file in enumerate([out_p1, out_p2]):
                ext = "/%s" % (i + 1)
                with file_transaction(out_file) as tx_out_file:
                    _merge_and_bgzip([xs[i] for xs in fastqs], tx_out_file, out_file, ext)
    if utils.file_exists(out_p1):
        return [out_p1, out_p2]
    else:
        assert utils.file_exists(out_s)
        return [out_s]

def _merge_and_bgzip(orig_files, out_file, base_file, ext=""):
    """Merge a group of gzipped input files into a final bgzipped output.

    Also handles providing unique names for each input file to avoid
    collisions on multi-region output. Handles renaming with awk magic from:
    https://www.biostars.org/p/68477/

    Removes orig_files after merging.
    """
    assert out_file.endswith(".gz")
    full_file = out_file.replace(".gz", "")
    run_file = "%s-merge.bash" % utils.splitext_plus(base_file)[0]

    cmds = ["set -e\n"]
    for i, fname in enumerate(orig_files):
        cmd = ("""zcat %s | awk '{print (NR%%4 == 1) ? "@%s_" ++i "%s" : $0}' >> %s\n"""
               % (fname, i, ext, full_file))
        cmds.append(cmd)
    cmds.append("bgzip %s\n" % full_file)

    with open(run_file, "w") as out_handle:
        out_handle.write("".join("".join(cmds)))
    do.run([do.find_bash(), run_file], "Rename, merge and bgzip CRAM fastq output")
    assert os.path.exists(out_file) and not _is_gzip_empty(out_file)
    for fname in orig_files:
        os.remove(fname)

def _cram_to_fastq_regions(regions, cram_file, dirs, data):
    """Convert CRAM files to fastq, potentially within sub regions.

    Returns multiple fastq files that can be merged back together.
    """
    base_name = utils.splitext_plus(os.path.basename(cram_file))[0]
    work_dir = utils.safe_makedir(os.path.join(dirs["work"], "align_prep",
                                               "%s-parts" % base_name))
    ref_file = tz.get_in(["reference", "fasta", "base"], data)
    resources = config_utils.get_resources("bamtofastq", data["config"])
    cores = tz.get_in(["config", "algorithm", "num_cores"], data, 1)
    max_mem = int(resources.get("memory", "1073741824")) * cores  # 1Gb/core default
    fnames = []
    is_paired = False
    for region in regions:
        rext = "-%s" % region.replace(":", "_").replace("-", "_") if region else "full"
        out_s, out_p1, out_p2 = [os.path.join(work_dir, "%s%s-%s.fq.gz" %
                                              (base_name, rext, fext))
                                 for fext in ["s1", "p1", "p2"]]
        if not utils.file_exists(out_p1):
            with file_transaction(out_s, out_p1, out_p2) as (tx_out_s, tx_out_p1, tx_out_p2):
                sortprefix = "%s-sort" % utils.splitext_plus(tx_out_s)[0]
                cmd = ("bamtofastq filename={cram_file} inputformat=cram T={sortprefix} "
                       "gz=1 collate=1 colsbs={max_mem} "
                       "F={tx_out_p1} F2={tx_out_p2} S={tx_out_s} O=/dev/null O2=/dev/null "
                       "reference={ref_file}")
                if region:
                    cmd += " ranges='{region}'"
                do.run(cmd.format(**locals()), "CRAM to fastq %s" % region if region else "")
        if is_paired or not _is_gzip_empty(out_p1):
            fnames.append((out_p1, out_p2))
            is_paired = True
        else:
            fnames.append((out_s,))
    return fnames

def _is_gzip_empty(fname):
    count = subprocess.check_output("zcat %s | head -1 | wc -l" % fname, shell=True,
                                    stderr=open("/dev/null", "w"))
    return int(count) < 1

def _bgzip_from_bam(bam_file, dirs, config, is_retry=False):
    """Create bgzipped fastq files from an input BAM file.
    """
    # tools
    bamtofastq = config_utils.get_program("bamtofastq", config)
    resources = config_utils.get_resources("bamtofastq", config)
    cores = config["algorithm"].get("num_cores", 1)
    max_mem = int(resources.get("memory", "1073741824")) * cores  # 1Gb/core default
    bgzip = tools.get_bgzip_cmd(config, is_retry)
    # files
    work_dir = utils.safe_makedir(os.path.join(dirs["work"], "align_prep"))
    out_file_1 = os.path.join(work_dir, "%s-1.fq.gz" % os.path.splitext(os.path.basename(bam_file))[0])
    if bam.is_paired(bam_file):
        out_file_2 = out_file_1.replace("-1.fq.gz", "-2.fq.gz")
    else:
        out_file_2 = None
    needs_retry = False
    if is_retry or not utils.file_exists(out_file_1):
        with file_transaction(out_file_1) as tx_out_file:
            for f in [tx_out_file, out_file_1, out_file_2]:
                if f and os.path.exists(f):
                    os.remove(f)
            fq1_bgzip_cmd = "%s -c /dev/stdin > %s" % (bgzip, tx_out_file)
            sortprefix = "%s-sort" % os.path.splitext(tx_out_file)[0]
            if bam.is_paired(bam_file):
                fq2_bgzip_cmd = "%s -c /dev/stdin > %s" % (bgzip, out_file_2)
                out_str = ("F=>({fq1_bgzip_cmd}) F2=>({fq2_bgzip_cmd}) S=/dev/null O=/dev/null "
                           "O2=/dev/null collate=1 colsbs={max_mem}")
            else:
                out_str = "S=>({fq1_bgzip_cmd})"
            cmd = "{bamtofastq} filename={bam_file} T={sortprefix} " + out_str
            try:
                do.run(cmd.format(**locals()), "BAM to bgzipped fastq",
                       checks=[do.file_reasonable_size(tx_out_file, bam_file)],
                       log_error=False)
            except subprocess.CalledProcessError, msg:
                if not is_retry and "deflate failed" in str(msg):
                    logger.info("bamtofastq deflate IO failure preparing %s. Retrying with single core."
                                % (bam_file))
                    needs_retry = True
                else:
                    logger.exception()
                    raise
    if needs_retry:
        return _bgzip_from_bam(bam_file, dirs, config, is_retry=True)
    else:
        return [x for x in [out_file_1, out_file_2] if x is not None]

@utils.map_wrap
@zeromq_aware_logging
def _grabix_index(data):
    in_file = data["bgzip_file"]
    config = data["config"]
    grabix = config_utils.get_program("grabix", config)
    gbi_file = in_file + ".gbi"
    if not utils.file_exists(gbi_file) or _is_partial_index(gbi_file):
        do.run([grabix, "index", in_file], "Index input with grabix: %s" % os.path.basename(in_file))
    return gbi_file

def _is_partial_index(gbi_file):
    """Check for truncated output since grabix doesn't write to a transactional directory.
    """
    with open(gbi_file) as in_handle:
        for i, _ in enumerate(in_handle):
            if i > 2:
                return False
    return True

@utils.map_wrap
@zeromq_aware_logging
def _bgzip_from_fastq(data):
    """Prepare a bgzipped file from a fastq input, potentially gzipped (or bgzipped already).
    """
    in_file = data["in_file"]
    config = data["config"]
    grabix = config_utils.get_program("grabix", config)
    needs_convert = config["algorithm"].get("quality_format", "").lower() == "illumina"
    if in_file.endswith(".gz"):
        needs_bgzip, needs_gunzip = _check_gzipped_input(in_file, grabix, needs_convert)
    else:
        needs_bgzip, needs_gunzip = True, False
    if needs_bgzip or needs_gunzip or needs_convert:
        out_file = _bgzip_file(in_file, data["dirs"], config, needs_bgzip, needs_gunzip,
                               needs_convert)
    else:
        out_file = in_file
    return [out_file]

def _bgzip_file(in_file, dirs, config, needs_bgzip, needs_gunzip, needs_convert):
    """Handle bgzip of input file, potentially gunzipping an existing file.
    """
    work_dir = utils.safe_makedir(os.path.join(dirs["work"], "align_prep"))
    out_file = os.path.join(work_dir, os.path.basename(in_file) +
                            (".gz" if not in_file.endswith(".gz") else ""))
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            assert needs_bgzip
            bgzip = tools.get_bgzip_cmd(config)
            if needs_convert:
                in_file = fastq_convert_pipe_cl(in_file, {"config": config})
            if needs_gunzip:
                gunzip_cmd = "gunzip -c {in_file} |".format(**locals())
                bgzip_in = "/dev/stdin"
            else:
                gunzip_cmd = ""
                bgzip_in = in_file
            do.run("{gunzip_cmd} {bgzip} -c {bgzip_in} > {tx_out_file}".format(**locals()),
                   "bgzip input file")
    return out_file

def _check_gzipped_input(in_file, grabix, needs_convert):
    """Determine if a gzipped input file is blocked gzip or standard.
    """
    is_bgzip = subprocess.check_output([grabix, "check", in_file])
    if is_bgzip.strip() == "yes" and not needs_convert:
        return False, False
    else:
        return True, True
