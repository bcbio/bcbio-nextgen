"""Prepare read inputs (fastq, gzipped fastq and BAM) for parallel NGS alignment.
"""
import collections
import copy
import glob
import os
import shutil
import subprocess

import six
import toolz as tz

from bcbio import bam, utils
from bcbio.bam import cram, fastq
from bcbio.cwl import cwlutils
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
    machine. Prepares a bgzip and grabix indexed file for retrieving sections
    of files.
    """
    from bcbio.pipeline import sample
    data = cwlutils.normalize_missing(data)
    aligner = tz.get_in(("config", "algorithm", "aligner"), data)
    # CRAM files must be converted to bgzipped fastq, unless not aligning.
    # Also need to prep and download remote files.
    if not ("files" in data and data["files"] and aligner and (_is_cram_input(data["files"]) or
                                                               objectstore.is_remote(data["files"][0]))):
        # skip indexing on samples without input files or not doing alignment
        if ("files" not in data or not data["files"] or data["files"][0] is None or not aligner):
            return [[data]]
    data["files_orig"] = data["files"]
    data["files"] = prep_fastq_inputs(data["files"], data)
    # preparation converts illumina into sanger format
    data["config"]["algorithm"]["quality_format"] = "standard"
    # Handle any necessary trimming
    data = utils.to_single_data(sample.trim_sample(data)[0])
    _prep_grabix_indexes(data["files"], data)
    data = _set_align_split_size(data)
    out = []
    if tz.get_in(["config", "algorithm", "align_split_size"], data):
        splits = _find_read_splits(data["files"][0], int(data["config"]["algorithm"]["align_split_size"]))
        for split in splits:
            cur_data = copy.deepcopy(data)
            cur_data["align_split"] = split
            out.append([cur_data])
    else:
        out.append([data])
    if "output_cwl_keys" in data:
        out = cwlutils.samples_to_records([utils.to_single_data(x) for x in out],
                                          ["files", "align_split", "config__algorithm__quality_format"])
    return out

def _set_align_split_size(data):
    """Set useful align_split_size, generating an estimate if it doesn't exist.

    We try to split on larger inputs and avoid too many pieces, aiming for size
    chunks of 5Gb or at most 50 maximum splits.

    The size estimate used in calculations is 20 million reads for ~5Gb.

    For UMI calculations we skip splitting since we're going to align and
    re-align after consensus.

    For CWL runs, we pick larger split sizes to avoid overhead of staging each chunk.
    """
    if cwlutils.is_cwl_run(data):
        target_size = 20  # Gb
        target_size_reads = 80  # million reads
    else:
        target_size = 5  # Gb
        target_size_reads = 20  # million reads
    max_splits = 100  # Avoid too many pieces, causing merge memory problems
    val = dd.get_align_split_size(data)
    umi_consensus = dd.get_umi_consensus(data)
    if val is None:
        if not umi_consensus:
            total_size = 0  # Gb
            # Use original files if we might have reduced the size of our prepped files
            input_files = data.get("files_orig", []) if dd.get_save_diskspace(data) else data.get("files", [])
            for fname in input_files:
                if os.path.exists(fname):
                    total_size += os.path.getsize(fname) / (1024.0 * 1024.0 * 1024.0)
            # Only set if we have files and are bigger than the target size
            if total_size > target_size:
                data["config"]["algorithm"]["align_split_size"] = \
                  int(1e6 * _pick_align_split_size(total_size, target_size,
                                                   target_size_reads, max_splits))
    elif val:
        assert not umi_consensus, "Cannot set align_split_size to %s with UMI conensus specified" % val
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
                assert _get_grabix_index(in_file), "Need grabix index for %s" % in_file
                out.append("<(grabix grab {in_file} {start} {end})".format(**locals()))
            else:
                out.append(None)
        return out

def _seqtk_fastq_prep_cl(data, in_file=None, read_num=0):
    """Provide a commandline for prep of fastq inputs with seqtk.

    Handles fast conversion of fastq quality scores and trimming.
    """
    needs_convert = dd.get_quality_format(data).lower() == "illumina"
    trim_ends = dd.get_trim_ends(data)
    seqtk = config_utils.get_program("seqtk", data["config"])
    if in_file:
        in_file = objectstore.cl_input(in_file)
    else:
        in_file = "/dev/stdin"
    cmd = ""
    if needs_convert:
        cmd += "{seqtk} seq -Q64 -V {in_file}".format(**locals())
    if trim_ends:
        left_trim, right_trim = trim_ends[0:2] if data.get("read_num", read_num) == 0 else trim_ends[2:4]
        if left_trim or right_trim:
            trim_infile = "/dev/stdin" if needs_convert else in_file
            pipe = " | " if needs_convert else ""
            cmd += "{pipe}{seqtk} trimfq -b {left_trim} -e {right_trim} {trim_infile}".format(**locals())
    return cmd

def fastq_convert_pipe_cl(in_file, data):
    """Create an anonymous pipe converting Illumina 1.3-1.7 to Sanger.

    Uses seqtk: https://github.com/lh3/seqt
    """
    cmd = _seqtk_fastq_prep_cl(data, in_file)
    if not cmd:
        cat_cmd = "zcat" if in_file.endswith(".gz") else "cat"
        cmd = cat_cmd + " " + in_file
    return "<(%s)" % cmd

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
    for mgroup in to_merge.values():
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
    if not tz.get_in(["config", "algorithm", "kraken"], data):
        # kraken requires fasta filenames from data['files'] as input.
        # We don't want to remove those files if kraken qc is required.
        _save_fastq_space(samples)
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
        else:
            hla_files = glob.glob(os.path.join(dd.get_work_dir(data), "align",
                                               dd.get_sample_name(data), "hla", "*.fq"))
            if hla_files:
                data["hla"]["fastq"] = hla_files
        out.append([data])
    return out

def _save_fastq_space(items):
    """Potentially save fastq space prior to merging, since alignments done.
    """
    to_cleanup = {}
    for data in (utils.to_single_data(x) for x in items):
        for fname in data.get("files", []):
            if os.path.realpath(fname).startswith(dd.get_work_dir(data)):
                to_cleanup[fname] = data["config"]
    for fname, config in to_cleanup.items():
        utils.save_diskspace(fname, "Cleanup prep files after alignment finished", config)

# ## determine file sections

def _get_grabix_index(in_file):
    gbi_file = in_file + ".gbi"
    if utils.file_exists(gbi_file):
        with open(gbi_file) as in_handle:
            header = next(in_handle)
            if header.find("Not grabix indexed") == -1:
                return gbi_file

def get_downsample_params(data):
    ds_mult = dd.get_maxcov_downsample(data)
    if ds_mult and ds_mult > 0:
        return {"min_coverage_for_downsampling": 10,
                "maxcov_downsample_multiplier": ds_mult}

def total_reads_from_grabix(in_file):
    """Retrieve total reads in a fastq file from grabix index.
    """
    gbi_file = _get_grabix_index(in_file)
    if gbi_file:
        with open(gbi_file) as in_handle:
            next(in_handle)  # throw away
            num_lines = int(next(in_handle).strip())
        assert num_lines % 4 == 0, "Expected lines to be multiple of 4"
        return num_lines // 4
    else:
        return 0

def _find_read_splits(in_file, split_size):
    """Determine sections of fastq files to process in splits.

    Assumes a 4 line order to input files (name, read, name, quality).
    grabix is 1-based inclusive, so return coordinates in that format.
    """
    num_lines = total_reads_from_grabix(in_file) * 4
    assert num_lines and num_lines > 0, "Did not find grabix index reads: %s %s" % (in_file, num_lines)
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

def _ready_gzip_fastq(in_files, data, require_bgzip=False):
    """Check if we have gzipped fastq and don't need format conversion or splitting.

    Avoid forcing bgzip if we don't need indexed files.
    """
    all_gzipped = all([not x or x.endswith(".gz") for x in in_files])
    if require_bgzip and all_gzipped:
        all_gzipped = all([not x or not _check_gzipped_input(x, data)[0] for x in in_files])
    needs_convert = dd.get_quality_format(data).lower() == "illumina"
    needs_trim = dd.get_trim_ends(data)
    do_splitting = dd.get_align_split_size(data) is not False
    return (all_gzipped and not needs_convert and not do_splitting and
            not objectstore.is_remote(in_files[0]) and not needs_trim and not get_downsample_params(data))

def prep_fastq_inputs(in_files, data):
    """Prepare bgzipped fastq inputs
    """
    if len(in_files) == 1 and _is_bam_input(in_files):
        out = _bgzip_from_bam(in_files[0], data["dirs"], data)
    elif len(in_files) == 1 and _is_cram_input(in_files):
        out = _bgzip_from_cram(in_files[0], data["dirs"], data)
    elif len(in_files) in [1, 2] and _ready_gzip_fastq(in_files, data):
        out = _symlink_in_files(in_files, data)
    else:
        if len(in_files) > 2:
            fpairs = fastq.combine_pairs(in_files)
            pair_types = set([len(xs) for xs in fpairs])
            assert len(pair_types) == 1
            fpairs.sort(key=lambda x: os.path.basename(x[0]))
            organized = [[xs[0] for xs in fpairs]]
            if len(fpairs[0]) > 1:
                organized.append([xs[1] for xs in fpairs])
            in_files = organized
        parallel = {"type": "local", "num_jobs": len(in_files),
                    "cores_per_job": max(1, data["config"]["algorithm"]["num_cores"] // len(in_files))}
        inputs = [{"in_file": x, "read_num": i, "dirs": data["dirs"], "config": data["config"],
                   "is_cwl": "cwl_keys" in data,
                   "rgnames": data["rgnames"]}
                  for i, x in enumerate(in_files) if x]
        out = run_multicore(_bgzip_from_fastq_parallel, [[d] for d in inputs], data["config"], parallel)
    return out

def _symlink_in_files(in_files, data):
    """Symlink (shared filesystem) or copy (CWL) inputs into align_prep directory.
    """
    work_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "align_prep"))
    out = []
    for in_file in in_files:
        out_file = os.path.join(work_dir, "%s_%s" % (dd.get_sample_name(data), os.path.basename(in_file)))
        out_file = _symlink_or_copy_grabix(in_file, out_file, data)
        out.append(out_file)
    return out

def _symlink_or_copy_grabix(in_file, out_file, data):
    """We cannot symlink in CWL, but may be able to use inputs or copy
    """
    if cwlutils.is_cwl_run(data):
        # Has grabix indexes, we're okay to go
        if utils.file_exists(in_file + ".gbi"):
            out_file = in_file
        else:
            utils.copy_plus(in_file, out_file)
    else:
        utils.symlink_plus(in_file, out_file)
    return out_file

def _prep_grabix_indexes(in_files, data):
    """Parallel preparation of grabix indexes for files.
    """
    # if we have gzipped but not bgzipped, add a fake index for CWL support
    # Also skips bgzip indexing if we don't need alignment splitting
    if _ready_gzip_fastq(in_files, data) and (not _ready_gzip_fastq(in_files, data, require_bgzip=True) or
                                              dd.get_align_split_size(data) is False):
        for in_file in in_files:
            if not utils.file_exists(in_file + ".gbi"):
                with file_transaction(data, in_file + ".gbi") as tx_gbi_file:
                    with open(tx_gbi_file, "w") as out_handle:
                        out_handle.write("Not grabix indexed; index added for compatibility.\n")
    else:
        items = [[{"bgzip_file": x, "config": copy.deepcopy(data["config"])}] for x in in_files if x]
        run_multicore(_grabix_index, items, data["config"])
    return data

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

def _bgzip_from_bam(bam_file, dirs, data, is_retry=False, output_infix=''):
    """Create bgzipped fastq files from an input BAM file.
    """
    # tools
    config = data["config"]
    bamtofastq = config_utils.get_program("bamtofastq", config)
    resources = config_utils.get_resources("bamtofastq", config)
    cores = config["algorithm"].get("num_cores", 1)
    max_mem = config_utils.convert_to_bytes(resources.get("memory", "1G")) * cores
    bgzip = tools.get_bgzip_cmd(config, is_retry)
    # files
    work_dir = utils.safe_makedir(os.path.join(dirs["work"], "align_prep"))
    out_file_1 = os.path.join(work_dir, "%s%s-1.fq.gz" % (os.path.splitext(os.path.basename(bam_file))[0], output_infix))
    out_file_2 = out_file_1.replace("-1.fq.gz", "-2.fq.gz")
    needs_retry = False
    if is_retry or not utils.file_exists(out_file_1):
        if not bam.is_paired(bam_file):
            out_file_2 = None
        with file_transaction(config, out_file_1) as tx_out_file:
            for f in [tx_out_file, out_file_1, out_file_2]:
                if f and os.path.exists(f):
                    os.remove(f)
            fq1_bgzip_cmd = "%s -c /dev/stdin > %s" % (bgzip, tx_out_file)
            prep_cmd = _seqtk_fastq_prep_cl(data, read_num=0)
            if prep_cmd:
                fq1_bgzip_cmd = prep_cmd + " | " + fq1_bgzip_cmd
            sortprefix = "%s-sort" % os.path.splitext(tx_out_file)[0]
            if bam.is_paired(bam_file):
                prep_cmd = _seqtk_fastq_prep_cl(data, read_num=1)
                fq2_bgzip_cmd = "%s -c /dev/stdin > %s" % (bgzip, out_file_2)
                if prep_cmd:
                    fq2_bgzip_cmd = prep_cmd + " | " + fq2_bgzip_cmd
                out_str = ("F=>({fq1_bgzip_cmd}) F2=>({fq2_bgzip_cmd}) S=/dev/null O=/dev/null "
                           "O2=/dev/null collate=1 colsbs={max_mem}")
            else:
                out_str = "S=>({fq1_bgzip_cmd})"
            bam_file = objectstore.cl_input(bam_file)
            extra_opts = " ".join([str(x) for x in resources.get("options", [])])
            cmd = "{bamtofastq} filename={bam_file} T={sortprefix} {extra_opts} " + out_str
            try:
                do.run(cmd.format(**locals()), "BAM to bgzipped fastq",
                       checks=[do.file_reasonable_size(tx_out_file, bam_file)],
                       log_error=False)
            except subprocess.CalledProcessError as msg:
                if not is_retry and "deflate failed" in str(msg):
                    logger.info("bamtofastq deflate IO failure preparing %s. Retrying with single core."
                                % (bam_file))
                    needs_retry = True
                else:
                    logger.exception()
                    raise
    if needs_retry:
        return _bgzip_from_bam(bam_file, dirs, data, is_retry=True)
    else:
        return [x for x in [out_file_1, out_file_2] if x is not None and utils.file_exists(x)]

@utils.map_wrap
@zeromq_aware_logging
def _grabix_index(data):
    """Create grabix index of bgzip input file.

    grabix does not allow specification of output file, so symlink the original
    file into a transactional directory.
    """
    in_file = data["bgzip_file"]
    config = data["config"]
    grabix = config_utils.get_program("grabix", config)
    gbi_file = _get_grabix_index(in_file)
    # We always build grabix input so we can use it for counting reads and doing downsampling
    if not gbi_file or _is_partial_index(gbi_file):
        if gbi_file:
            utils.remove_safe(gbi_file)
        else:
            gbi_file = in_file + ".gbi"
        with file_transaction(data, gbi_file) as tx_gbi_file:
            tx_in_file = os.path.splitext(tx_gbi_file)[0]
            utils.symlink_plus(in_file, tx_in_file)
            do.run([grabix, "index", tx_in_file], "Index input with grabix: %s" % os.path.basename(in_file))
    assert utils.file_exists(gbi_file)
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
    if isinstance(in_file, (list, tuple)):
        in_file = in_file[0]
    needs_convert = dd.get_quality_format(data).lower() == "illumina"
    # special case, empty files that have been cleaned
    if not objectstore.is_remote(in_file) and os.path.getsize(in_file) == 0:
        needs_bgzip, needs_gunzip = False, False
    elif in_file.endswith(".gz") and not objectstore.is_remote(in_file):
        if needs_convert or dd.get_trim_ends(data):
            needs_bgzip, needs_gunzip = True, True
        else:
            needs_bgzip, needs_gunzip = _check_gzipped_input(in_file, data)
    elif in_file.endswith(".bz2"):
        needs_bgzip, needs_gunzip = True, True
    elif objectstore.is_remote(in_file) and not tz.get_in(["config", "algorithm", "align_split_size"], data):
        needs_bgzip, needs_gunzip = False, False
    else:
        needs_bgzip, needs_gunzip = True, False
    work_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "align_prep"))
    if (needs_bgzip or needs_gunzip or needs_convert or dd.get_trim_ends(data) or
          objectstore.is_remote(in_file) or
          (isinstance(data["in_file"], (tuple, list)) and len(data["in_file"]) > 1)):
        out_file = _bgzip_file(data["in_file"], data["config"], work_dir,
                               needs_bgzip, needs_gunzip, needs_convert, data)
    else:
        out_file = os.path.join(work_dir, "%s_%s" % (dd.get_sample_name(data), os.path.basename(in_file)))
        out_file = _symlink_or_copy_grabix(in_file, out_file, data)
    return out_file

def _bgzip_file(finput, config, work_dir, needs_bgzip, needs_gunzip, needs_convert, data):
    """Handle bgzip of input file, potentially gunzipping an existing file.

    Handles cases where finput might be multiple files and need to be concatenated.
    """
    if isinstance(finput, six.string_types):
        in_file = finput
    else:
        assert not needs_convert, "Do not yet handle quality conversion with multiple inputs"
        return _bgzip_multiple_files(finput, work_dir, data)
    out_file = os.path.join(work_dir, os.path.basename(in_file).replace(".bz2", "") +
                            (".gz" if not in_file.endswith(".gz") else ""))
    if not utils.file_exists(out_file):
        with file_transaction(config, out_file) as tx_out_file:
            bgzip = tools.get_bgzip_cmd(config)
            is_remote = objectstore.is_remote(in_file)
            in_file = objectstore.cl_input(in_file, unpack=needs_gunzip or needs_convert or
                                           needs_bgzip or dd.get_trim_ends(data))
            if needs_convert or dd.get_trim_ends(data):
                in_file = fastq_convert_pipe_cl(in_file, data)
            if needs_gunzip and not (needs_convert or dd.get_trim_ends(data)):
                if in_file.endswith(".bz2"):
                    gunzip_cmd = "bunzip2 -c {in_file} |".format(**locals())
                else:
                    gunzip_cmd = "gunzip -c {in_file} |".format(**locals())
                bgzip_in = "/dev/stdin"
            else:
                gunzip_cmd = ""
                bgzip_in = in_file
            if needs_bgzip:
                do.run("{gunzip_cmd} {bgzip} -c {bgzip_in} > {tx_out_file}".format(**locals()),
                       "bgzip input file")
            elif is_remote:
                bgzip = "| bgzip -c" if (needs_convert or dd.get_trim_ends(data)) else ""
                do.run("cat {in_file} {bgzip} > {tx_out_file}".format(**locals()), "Get remote input")
            else:
                raise ValueError("Unexpected inputs: %s %s %s %s" % (in_file, needs_bgzip,
                                                                     needs_gunzip, needs_convert))
    return out_file

def _bgzip_multiple_files(in_files, work_dir, data):
    out_file = os.path.join(work_dir, "%s-combined-%s" % (dd.get_sample_name(data),
                                                          os.path.basename(in_files[0]).replace(".bz2", "") +
                                                          (".gz" if not in_files[0].endswith(".gz") else "")))
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            if in_files[0].endswith(".bz2"):
                gunzip_cmd = "bunzip2 -c"
            elif in_files[0].endswith(".gz"):
                gunzip_cmd = "gunzip -c"
            else:
                gunzip_cmd = "cat"
            cmd = "%s %s | bgzip -c > %s" % (gunzip_cmd, " ".join(in_files), tx_out_file)
            do.run(cmd, "Combine and bgzip multiple input files: %s" % dd.get_sample_name(data))
    return out_file

def _check_gzipped_input(in_file, data):
    """Determine if a gzipped input file is blocked gzip or standard.
    """
    grabix = config_utils.get_program("grabix", data["config"])
    is_bgzip = subprocess.check_output([grabix, "check", in_file])
    if is_bgzip.strip() == "yes":
        return False, False
    else:
        return True, True
