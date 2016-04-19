"""Prepare read inputs (fastq, gzipped fastq and BAM) for parallel NGS alignment.
"""
import collections
import copy
import glob
import os
import shutil
import subprocess

import toolz as tz

from bcbio import bam, utils
from bcbio.bam import cram
from bcbio.log import logger
from bcbio.distributed import objectstore
from bcbio.distributed.multi import run_multicore, zeromq_aware_logging
from bcbio.distributed.transaction import file_transaction
from bcbio.ngsalign import rtg
from bcbio.pipeline import config_utils, tools
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do

def create_inputs(data):
    """Index input reads and prepare groups of reads to process concurrently.

    Allows parallelization of alignment beyond processors available on a single
    machine. Prepares a rtg SDF format file with build in indexes for retrieving
    sections of files.

    Retains back compatibility with bgzip/grabix approach.
    """
    aligner = tz.get_in(("config", "algorithm", "aligner"), data)
    # CRAM files must be converted to bgzipped fastq, unless not aligning.
    # Also need to prep and download remote files.
    if not ("files" in data and data["files"] and aligner and (_is_cram_input(data["files"]) or
                                                               objectstore.is_remote(data["files"][0]))):
        # skip indexing on samples without input files or not doing alignment
        if ("files" not in data or not data["files"] or data["files"][0] is None or not aligner):
            return [[data]]
    approach = "grabix" if _has_grabix_indices(data) else dd.get_align_prep_method(data)
    if approach == "rtg":
        data["files"] = [rtg.to_sdf(data["files"], data)]
    else:
        data["files"] = _prep_grabix_indexes(data["files"], data["dirs"], data)
    # preparation converts illumina into sanger format
    data["config"]["algorithm"]["quality_format"] = "standard"
    data = _set_align_split_size(data)
    if tz.get_in(["config", "algorithm", "align_split_size"], data):
        out = []
        if approach == "rtg":
            splits = rtg.calculate_splits(data["files"][0], data["config"]["algorithm"]["align_split_size"])
        else:
            splits = _find_read_splits(data["files"][0], data["config"]["algorithm"]["align_split_size"])
        for split in splits:
            cur_data = copy.deepcopy(data)
            cur_data["align_split"] = split
            out.append([cur_data])
        return out
    else:
        return [[data]]

def _set_align_split_size(data):
    """Set useful align_split_size, generating an estimate if it doesn't exist.

    We try to split on larger inputs and avoid too many pieces, aiming for size
    chunks of 5Gb or at most 50 maximum splits.

    The size estimate used in calculations is 20 million reads for ~5Gb.
    """
    target_size = 5  # Gb
    target_size_reads = 20  # million reads
    max_splits = 100  # Avoid too many pieces, causing merge memory problems
    val = tz.get_in(["config", "algorithm", "align_split_size"], data)
    if val is None:
        total_size = 0  # Gb
        for fname in data.get("files", []):
            if os.path.exists(fname):
                total_size += os.path.getsize(fname) / (1024.0 * 1024.0 * 1024.0)
        # Only set if we have files and are bigger than the target size
        if total_size > target_size:
            data["config"]["algorithm"]["align_split_size"] = \
              int(1e6 * _pick_align_split_size(total_size, target_size,
                                               target_size_reads, max_splits))
    return data

def _pick_align_split_size(total_size, target_size, target_size_reads, max_splits):
    """Do the work of picking an alignment split size for the given criteria.
    """
    # Too many pieces, increase our target size to get max_splits pieces
    if total_size // target_size > max_splits:
        piece_size = total_size // max_splits
        return int(piece_size * target_size_reads / target_size)
    else:
        return int(target_size_reads)

def _has_grabix_indices(data):
    """Back compatibility with existing runs, look for grabix indexes.
    """
    work_dir = (os.path.join(data["dirs"]["work"], "align_prep"))
    return len(glob.glob(os.path.join(work_dir, "*.gbi"))) > 0

def split_namedpipe_cls(pair1_file, pair2_file, data):
    """Create a commandline suitable for use as a named pipe with reads in a given region.
    """
    if "align_split" in data:
        start, end = [int(x) for x in data["align_split"].split("-")]
    else:
        start, end = None, None
    if pair1_file.endswith(".sdf"):
        assert not pair2_file, pair2_file
        return rtg.to_fastq_apipe_cl(pair1_file, start, end)
    else:
        out = []
        for in_file in pair1_file, pair2_file:
            if in_file:
                assert os.path.exists(in_file + ".gbi"), "Need grabix index for %s" % in_file
                out.append("<(grabix grab {in_file} {start} {end})".format(**locals()))
            else:
                out.append(None)
        return out

def fastq_convert_pipe_cl(in_file, data):
    """Create an anonymous pipe converting Illumina 1.3-1.7 to Sanger.

    Uses seqtk: https://github.com/lh3/seqt
    """
    seqtk = config_utils.get_program("seqtk", data["config"])
    in_file = objectstore.cl_input(in_file)
    return "<({seqtk} seq -Q64 -V {in_file})".format(**locals())

# ## configuration

def parallel_multiplier(items):
    """Determine if we will be parallelizing items during processing.
    """
    multiplier = 1
    for data in (x[0] for x in items):
        if (tz.get_in(["config", "algorithm", "align_split_size"], data) is not False and
              tz.get_in(["algorithm", "align_split_size"], data) is not False):
            multiplier += 50
    return multiplier

# ## merge

def setup_combine(final_file, data):
    """Setup the data and outputs to allow merging data back together.
    """
    if "align_split" not in data:
        return final_file, data
    align_dir = os.path.dirname(final_file)
    base, ext = os.path.splitext(os.path.basename(final_file))
    start, end = [int(x) for x in data["align_split"].split("-")]
    out_file = os.path.join(utils.safe_makedir(os.path.join(align_dir, "split")),
                            "%s-%s_%s%s" % (base, start, end, ext))
    data["combine"] = {"work_bam": {"out": final_file, "extras": []}}
    return out_file, data

def merge_split_alignments(samples, run_parallel):
    """Manage merging split alignments back into a final working BAM file.

    Perform de-duplication on the final merged file.
    """
    ready = []
    file_key = "work_bam"
    to_merge = collections.defaultdict(list)
    for data in (xs[0] for xs in samples):
        if data.get("combine"):
            out_key = tz.get_in(["combine", file_key, "out"], data)
            if not out_key:
                out_key = data["rgnames"]["lane"]
            to_merge[out_key].append(data)
        else:
            ready.append([data])
    ready_merge = []
    hla_merges = []
    for mgroup in to_merge.itervalues():
        cur_data = mgroup[0]
        del cur_data["align_split"]
        for x in mgroup[1:]:
            cur_data["combine"][file_key]["extras"].append(x[file_key])
        ready_merge.append([cur_data])
        cur_hla = None
        for d in mgroup:
            hla_files = tz.get_in(["hla", "fastq"], d)
            if hla_files:
                if not cur_hla:
                    cur_hla = {"rgnames": {"sample": dd.get_sample_name(cur_data)},
                               "config": cur_data["config"], "dirs": cur_data["dirs"],
                               "hla": {"fastq": []}}
                cur_hla["hla"]["fastq"].append(hla_files)
        if cur_hla:
            hla_merges.append([cur_hla])
    merged = run_parallel("delayed_bam_merge", ready_merge)
    hla_merge_raw = run_parallel("merge_split_alignments", hla_merges)
    hla_merges = {}
    for hla_merge in [x[0] for x in hla_merge_raw]:
        hla_merges[dd.get_sample_name(hla_merge)] = tz.get_in(["hla", "fastq"], hla_merge)
    # Add stable 'align_bam' target to use for retrieving raw alignment
    out = []
    for data in [x[0] for x in merged + ready]:
        if data.get("work_bam"):
            data["align_bam"] = data["work_bam"]
        if dd.get_sample_name(data) in hla_merges:
            data["hla"]["fastq"] = hla_merges[dd.get_sample_name(data)]
        out.append([data])
    return out

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
        chunks.append((last, min(new, num_lines)))
        last = new + 1
    return ["%s-%s" % (s, e) for s, e in chunks]

# ## bgzip and grabix

def _is_bam_input(in_files):
    return in_files and in_files[0].endswith(".bam") and (len(in_files) == 1 or in_files[1] is None)

def _is_cram_input(in_files):
    return in_files and in_files[0].endswith(".cram") and (len(in_files) == 1 or in_files[1] is None)

def _ready_gzip_fastq(in_files, data):
    """Check if we have gzipped fastq and don't need format conversion or splitting.
    """
    all_gzipped = all([not x or x.endswith(".gz") for x in in_files])
    needs_convert = tz.get_in(["config", "algorithm", "quality_format"], data, "").lower() == "illumina"
    do_splitting = tz.get_in(["config", "algorithm", "align_split_size"], data) is not False
    return all_gzipped and not needs_convert and not do_splitting and not objectstore.is_remote(in_files[0])

def _prep_grabix_indexes(in_files, dirs, data):
    if _is_bam_input(in_files):
        out = _bgzip_from_bam(in_files[0], dirs, data["config"])
    elif _is_cram_input(in_files):
        out = _bgzip_from_cram(in_files[0], dirs, data)
    elif _ready_gzip_fastq(in_files, data):
        out = in_files
    else:
        inputs = [{"in_file": x, "dirs": dirs, "config": data["config"], "rgnames": data["rgnames"]}
                  for x in in_files if x]
        if "pbgzip" not in dd.get_tools_off(data):
            out = [_bgzip_from_fastq(d) for d in inputs]
        else:
            out = run_multicore(_bgzip_from_fastq_parallel, [[d] for d in inputs], data["config"])
    items = [[{"bgzip_file": x, "config": copy.deepcopy(data["config"])}] for x in out if x]
    run_multicore(_grabix_index, items, data["config"])
    return out

def _bgzip_from_cram(cram_file, dirs, data):
    """Create bgzipped fastq files from an input CRAM file in regions of interest.

    Returns a list with a single file, for single end CRAM files, or two
    files for paired end input.
    """
    import pybedtools
    region_file = (tz.get_in(["config", "algorithm", "variant_regions"], data)
                   if tz.get_in(["config", "algorithm", "coverage_interval"], data)
                     in ["regional", "exome", "amplicon"]
                   else None)
    if region_file:
        regions = ["%s:%s-%s" % tuple(r[:3]) for r in pybedtools.BedTool(region_file)]
    else:
        regions = [None]
    work_dir = utils.safe_makedir(os.path.join(dirs["work"], "align_prep"))
    out_s, out_p1, out_p2 = [os.path.join(work_dir, "%s-%s.fq.gz" %
                                          (utils.splitext_plus(os.path.basename(cram_file))[0], fext))
                             for fext in ["s1", "p1", "p2"]]
    if (not utils.file_exists(out_s) and
          (not utils.file_exists(out_p1) or not utils.file_exists(out_p2))):
        cram.index(cram_file, data["config"])
        fastqs, part_dir = _cram_to_fastq_regions(regions, cram_file, dirs, data)
        if len(fastqs[0]) == 1:
            with file_transaction(data, out_s) as tx_out_file:
                _merge_and_bgzip([xs[0] for xs in fastqs], tx_out_file, out_s)
        else:
            for i, out_file in enumerate([out_p1, out_p2]):
                if not utils.file_exists(out_file):
                    ext = "/%s" % (i + 1)
                    with file_transaction(data, out_file) as tx_out_file:
                        _merge_and_bgzip([xs[i] for xs in fastqs], tx_out_file, out_file, ext)
        shutil.rmtree(part_dir)
    if utils.file_exists(out_p1):
        return [out_p1, out_p2]
    else:
        assert utils.file_exists(out_s)
        return [out_s]

def _bgzip_from_cram_sambamba(cram_file, dirs, data):
    """Use sambamba to extract from CRAM via regions.
    """
    raise NotImplementedError("sambamba doesn't yet support retrieval from CRAM by BED file")
    region_file = (tz.get_in(["config", "algorithm", "variant_regions"], data)
                   if tz.get_in(["config", "algorithm", "coverage_interval"], data) in ["regional", "exome"]
                   else None)
    base_name = utils.splitext_plus(os.path.basename(cram_file))[0]
    work_dir = utils.safe_makedir(os.path.join(dirs["work"], "align_prep",
                                               "%s-parts" % base_name))
    f1, f2, o1, o2, si = [os.path.join(work_dir, "%s.fq" % x) for x in ["match1", "match2", "unmatch1", "unmatch2",
                                                                        "single"]]
    ref_file = dd.get_ref_file(data)
    region = "-L %s" % region_file if region_file else ""
    cmd = ("sambamba view -f bam -l 0 -C {cram_file} -T {ref_file} {region} | "
           "bamtofastq F={f1} F2={f2} S={si} O={o1} O2={o2}")
    do.run(cmd.format(**locals()), "Convert CRAM to fastq in regions")

def _merge_and_bgzip(orig_files, out_file, base_file, ext=""):
    """Merge a group of gzipped input files into a final bgzipped output.

    Also handles providing unique names for each input file to avoid
    collisions on multi-region output. Handles renaming with awk magic from:
    https://www.biostars.org/p/68477/
    """
    assert out_file.endswith(".gz")
    full_file = out_file.replace(".gz", "")
    run_file = "%s-merge.bash" % utils.splitext_plus(base_file)[0]

    cmds = ["set -e\n"]
    for i, fname in enumerate(orig_files):
        cmd = ("""zcat %s | awk '{print (NR%%4 == 1) ? "@%s_" ++i "%s" : $0}' >> %s\n"""
               % (fname, i, ext, full_file))
        cmds.append(cmd)
    cmds.append("bgzip -f %s\n" % full_file)

    with open(run_file, "w") as out_handle:
        out_handle.write("".join("".join(cmds)))
    do.run([do.find_bash(), run_file], "Rename, merge and bgzip CRAM fastq output")
    assert os.path.exists(out_file) and not _is_gzip_empty(out_file)

def _cram_to_fastq_regions(regions, cram_file, dirs, data):
    """Convert CRAM files to fastq, potentially within sub regions.

    Returns multiple fastq files that can be merged back together.
    """
    base_name = utils.splitext_plus(os.path.basename(cram_file))[0]
    work_dir = utils.safe_makedir(os.path.join(dirs["work"], "align_prep",
                                               "%s-parts" % base_name))
    fnames = run_multicore(_cram_to_fastq_region,
                           [(cram_file, work_dir, base_name, region, data) for region in regions],
                           data["config"])
    # check if we have paired or single end data
    if any(not _is_gzip_empty(p1) for p1, p2, s in fnames):
        out = [[p1, p2] for p1, p2, s in fnames]
    else:
        out = [[s] for p1, p2, s in fnames]
    return out, work_dir

@utils.map_wrap
@zeromq_aware_logging
def _cram_to_fastq_region(cram_file, work_dir, base_name, region, data):
    """Convert CRAM to fastq in a specified region.
    """
    ref_file = tz.get_in(["reference", "fasta", "base"], data)
    resources = config_utils.get_resources("bamtofastq", data["config"])
    cores = tz.get_in(["config", "algorithm", "num_cores"], data, 1)
    max_mem = config_utils.convert_to_bytes(resources.get("memory", "1G")) * cores
    rext = "-%s" % region.replace(":", "_").replace("-", "_") if region else "full"
    out_s, out_p1, out_p2, out_o1, out_o2 = [os.path.join(work_dir, "%s%s-%s.fq.gz" %
                                                          (base_name, rext, fext))
                                             for fext in ["s1", "p1", "p2", "o1", "o2"]]
    if not utils.file_exists(out_p1):
        with file_transaction(data, out_s, out_p1, out_p2, out_o1, out_o2) as \
             (tx_out_s, tx_out_p1, tx_out_p2, tx_out_o1, tx_out_o2):
            cram_file = objectstore.cl_input(cram_file)
            sortprefix = "%s-sort" % utils.splitext_plus(tx_out_s)[0]
            cmd = ("bamtofastq filename={cram_file} inputformat=cram T={sortprefix} "
                   "gz=1 collate=1 colsbs={max_mem} exclude=SECONDARY,SUPPLEMENTARY "
                   "F={tx_out_p1} F2={tx_out_p2} S={tx_out_s} O={tx_out_o1} O2={tx_out_o2} "
                   "reference={ref_file}")
            if region:
                cmd += " ranges='{region}'"
            do.run(cmd.format(**locals()), "CRAM to fastq %s" % region if region else "")
    return [[out_p1, out_p2, out_s]]

def _is_gzip_empty(fname):
    count = subprocess.check_output("zcat %s | head -1 | wc -l" % fname, shell=True,
                                    stderr=open("/dev/null", "w"))
    return int(count) < 1

def _bgzip_from_bam(bam_file, dirs, config, is_retry=False, output_infix=''):
    """Create bgzipped fastq files from an input BAM file.
    """
    # tools
    bamtofastq = config_utils.get_program("bamtofastq", config)
    resources = config_utils.get_resources("bamtofastq", config)
    cores = config["algorithm"].get("num_cores", 1)
    max_mem = config_utils.convert_to_bytes(resources.get("memory", "1G")) * cores
    bgzip = tools.get_bgzip_cmd(config, is_retry)
    # files
    work_dir = utils.safe_makedir(os.path.join(dirs["work"], "align_prep"))
    out_file_1 = os.path.join(work_dir, "%s%s-1.fq.gz" % (os.path.splitext(os.path.basename(bam_file))[0], output_infix))
    if bam.is_paired(bam_file):
        out_file_2 = out_file_1.replace("-1.fq.gz", "-2.fq.gz")
    else:
        out_file_2 = None
    needs_retry = False
    if is_retry or not utils.file_exists(out_file_1):
        with file_transaction(config, out_file_1) as tx_out_file:
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
            bam_file = objectstore.cl_input(bam_file)
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
    if tz.get_in(["algorithm", "align_split_size"], config) is not False:
        if not utils.file_exists(gbi_file) or _is_partial_index(gbi_file):
            do.run([grabix, "index", in_file], "Index input with grabix: %s" % os.path.basename(in_file))
    return [gbi_file]

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
def _bgzip_from_fastq_parallel(data):
    return [_bgzip_from_fastq(data)]

def _bgzip_from_fastq(data):
    """Prepare a bgzipped file from a fastq input, potentially gzipped (or bgzipped already).
    """
    in_file = data["in_file"]
    config = data["config"]
    grabix = config_utils.get_program("grabix", config)
    needs_convert = config["algorithm"].get("quality_format", "").lower() == "illumina"
    if in_file.endswith(".gz") and not objectstore.is_remote(in_file):
        needs_bgzip, needs_gunzip = _check_gzipped_input(in_file, grabix, needs_convert)
    elif objectstore.is_remote(in_file) and not tz.get_in(["algorithm", "align_split_size"], config):
        needs_bgzip, needs_gunzip = False, False
    else:
        needs_bgzip, needs_gunzip = True, False
    work_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "align_prep"))
    if needs_bgzip or needs_gunzip or needs_convert or objectstore.is_remote(in_file):
        out_file = _bgzip_file(in_file, config, work_dir,
                               needs_bgzip, needs_gunzip, needs_convert)
    else:
        out_file = os.path.join(work_dir, "%s_%s" % (dd.get_sample_name(data), os.path.basename(in_file)))
        utils.symlink_plus(in_file, out_file)
    return out_file

def _bgzip_file(in_file, config, work_dir, needs_bgzip, needs_gunzip, needs_convert):
    """Handle bgzip of input file, potentially gunzipping an existing file.
    """
    out_file = os.path.join(work_dir, os.path.basename(in_file) +
                            (".gz" if not in_file.endswith(".gz") else ""))
    if not utils.file_exists(out_file):
        with file_transaction(config, out_file) as tx_out_file:
            bgzip = tools.get_bgzip_cmd(config)
            is_remote = objectstore.is_remote(in_file)
            in_file = objectstore.cl_input(in_file, unpack=needs_gunzip or needs_convert or needs_bgzip)
            if needs_convert:
                in_file = fastq_convert_pipe_cl(in_file, {"config": config})
            if needs_gunzip and not needs_convert:
                gunzip_cmd = "gunzip -c {in_file} |".format(**locals())
                bgzip_in = "/dev/stdin"
            else:
                gunzip_cmd = ""
                bgzip_in = in_file
            if needs_bgzip:
                do.run("{gunzip_cmd} {bgzip} -c {bgzip_in} > {tx_out_file}".format(**locals()),
                       "bgzip input file")
            elif is_remote:
                bgzip = "| bgzip -c" if needs_convert else ""
                do.run("cat {in_file} {bgzip} > {tx_out_file}".format(**locals()), "Get remote input")
            else:
                raise ValueError("Unexpected inputs: %s %s %s %s" % (in_file, needs_bgzip,
                                                                     needs_gunzip, needs_convert))
    return out_file

def _check_gzipped_input(in_file, grabix, needs_convert):
    """Determine if a gzipped input file is blocked gzip or standard.

    """
    # grabix is not parsing bgzip output from FASTQ correctly and will hang
    # indexing indefinitely. Until we can identify the issue, we re-bgzip everyting
    return True, True
    # is_bgzip = subprocess.check_output([grabix, "check", in_file])
    # if is_bgzip.strip() == "yes" and not needs_convert:
    #     return False, False
    # else:
    #     return True, True
