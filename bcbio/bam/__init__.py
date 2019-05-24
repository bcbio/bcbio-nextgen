"""Functionality to query and extract information from aligned BAM files.
"""
from __future__ import print_function
import collections
import os
import signal
import subprocess
import numpy

import pybedtools
import pysam
import toolz as tz
from six.moves import zip_longest

from bcbio import broad, utils
from bcbio.bam import ref
from bcbio.distributed import objectstore
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import config_utils
import bcbio.pipeline.datadict as dd
from bcbio.provenance import do

def is_empty(bam_file):
    """Determine if a BAM file is empty
    """
    bam_file = objectstore.cl_input(bam_file)
    cmd = ("set -o pipefail; "
           "samtools view {bam_file} | head -1 | wc -l")
    p = subprocess.Popen(cmd.format(**locals()), shell=True,
                         executable=do.find_bash(),
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
    stdout, stderr = p.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()
    if ((p.returncode == 0 or p.returncode == 141) and
         (stderr == "" or (stderr.startswith("gof3r") and stderr.endswith("broken pipe")))):
        return int(stdout) == 0
    else:
        raise ValueError("Failed to check empty status of BAM file: %s" % str(stderr))

def is_paired(bam_file):
    """Determine if a BAM file has paired reads.

    Works around issues with head closing the samtools pipe using signal trick from:
    http://stackoverflow.com/a/12451083/252589
    """
    bam_file = objectstore.cl_input(bam_file)
    cmd = ("set -o pipefail; "
           "samtools view -h {bam_file} | head -300000 | "
           "samtools view -S -f 1 /dev/stdin  | head -1 | wc -l")
    p = subprocess.Popen(cmd.format(**locals()), shell=True,
                         executable=do.find_bash(),
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
    stdout, stderr = p.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()
    stderr = stderr.strip()
    if ((p.returncode == 0 or p.returncode == 141) and
         (stderr == "" or (stderr.startswith("gof3r") and stderr.endswith("broken pipe")))):
        return int(stdout) > 0
    else:
        raise ValueError("Failed to check paired status of BAM file: %s" % str(stderr))

def fake_index(in_bam, data):
    """Create a fake index file for namesorted BAMs. bais require by CWL for consistency.
    """
    index_file = "%s.bai" % in_bam
    if not utils.file_exists(index_file):
        with file_transaction(data, index_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                out_handle.write("name sorted -- no index")
    return index_file

def index(in_bam, config, check_timestamp=True):
    """Index a BAM file, skipping if index present.

    Centralizes BAM indexing providing ability to switch indexing approaches.
    """
    assert is_bam(in_bam), "%s in not a BAM file" % in_bam
    index_file = "%s.bai" % in_bam
    alt_index_file = "%s.bai" % os.path.splitext(in_bam)[0]
    if check_timestamp:
        bai_exists = utils.file_uptodate(index_file, in_bam) or utils.file_uptodate(alt_index_file, in_bam)
    else:
        bai_exists = utils.file_exists(index_file) or utils.file_exists(alt_index_file)
    if not bai_exists:
        # Remove old index files and re-run to prevent linking into tx directory
        for fname in [index_file, alt_index_file]:
            utils.remove_safe(fname)
        samtools = config_utils.get_program("samtools", config)
        num_cores = config["algorithm"].get("num_cores", 1)
        with file_transaction(config, index_file) as tx_index_file:
            cmd = "{samtools} index -@ {num_cores} {in_bam} {tx_index_file}"
            do.run(cmd.format(**locals()), "Index BAM file: %s" % os.path.basename(in_bam))
    return index_file if utils.file_exists(index_file) else alt_index_file

def remove(in_bam):
    """
    remove bam file and the index if exists
    """
    if utils.file_exists(in_bam):
        utils.remove_safe(in_bam)
    if utils.file_exists(in_bam + ".bai"):
        utils.remove_safe(in_bam + ".bai")

def idxstats(in_bam, data):
    """Return BAM index stats for the given file, using samtools idxstats.
    """
    index(in_bam, data["config"], check_timestamp=False)
    AlignInfo = collections.namedtuple("AlignInfo", ["contig", "length", "aligned", "unaligned"])
    samtools = config_utils.get_program("samtools", data["config"])
    idxstats_out = subprocess.check_output([samtools, "idxstats", in_bam]).decode()
    out = []
    for line in idxstats_out.split("\n"):
        if line.strip():
            contig, length, aligned, unaligned = line.split("\t")
            out.append(AlignInfo(contig, int(length), int(aligned), int(unaligned)))
    return out

def fai_from_bam(ref_file, bam_file, out_file, data):
    """Create a fai index with only contigs in the input BAM file.
    """
    contigs = set([x.contig for x in idxstats(bam_file, data)])
    if not utils.file_uptodate(out_file, bam_file):
        with open(ref.fasta_idx(ref_file, data["config"])) as in_handle:
            with file_transaction(data, out_file) as tx_out_file:
                with open(tx_out_file, "w") as out_handle:
                    for line in (l for l in in_handle if l.strip()):
                        if line.split()[0] in contigs:
                            out_handle.write(line)
    return out_file

def ref_file_from_bam(bam_file, data):
    """Subset a fasta input file to only a fraction of input contigs.
    """
    new_ref = os.path.join(utils.safe_makedir(os.path.join(dd.get_work_dir(data), "inputs", "ref")),
                           "%s-subset.fa" % dd.get_genome_build(data))
    if not utils.file_exists(new_ref):
        with file_transaction(data, new_ref) as tx_out_file:
            contig_file = "%s-contigs.txt" % utils.splitext_plus(new_ref)[0]
            with open(contig_file, "w") as out_handle:
                for contig in [x.contig for x in idxstats(bam_file, data) if x.contig != "*"]:
                    out_handle.write("%s\n" % contig)
            cmd = "seqtk subseq -l 100 %s %s > %s" % (dd.get_ref_file(data), contig_file, tx_out_file)
            do.run(cmd, "Subset %s to BAM file contigs" % dd.get_genome_build(data))
    ref.fasta_idx(new_ref, data["config"])
    runner = broad.runner_from_path("picard", data["config"])
    runner.run_fn("picard_index_ref", new_ref)
    return {"base": new_ref}

def get_downsample_pct(in_bam, target_counts, data):
    """Retrieve percentage of file to downsample to get to target counts.

    Avoids minimal downsample which is not especially useful for
    improving QC times; 90& or more of reads.
    """
    total = sum(x.aligned for x in idxstats(in_bam, data))
    with pysam.Samfile(in_bam, "rb") as work_bam:
        n_rgs = max(1, len(work_bam.header.get("RG", [])))
    rg_target = n_rgs * target_counts
    if total > rg_target:
        pct = float(rg_target) / float(total)
        if pct < 0.9:
            return pct

def get_aligned_reads(in_bam, data):
    index(in_bam, data["config"], check_timestamp=False)
    bam_stats = idxstats(in_bam, data)
    align = sum(x.aligned for x in bam_stats)
    unaligned = sum(x.unaligned for x in bam_stats)
    total = float(align + unaligned)
    return 1.0 * align / total

def downsample(in_bam, data, target_counts, work_dir=None):
    """Downsample a BAM file to the specified number of target counts.
    """
    index(in_bam, data["config"], check_timestamp=False)
    ds_pct = get_downsample_pct(in_bam, target_counts, data)
    if ds_pct:
        out_file = "%s-downsample%s" % os.path.splitext(in_bam)
        if work_dir:
            out_file = os.path.join(work_dir, os.path.basename(out_file))
        if not utils.file_exists(out_file):
            with file_transaction(data, out_file) as tx_out_file:
                samtools = config_utils.get_program("samtools", data["config"])
                num_cores = dd.get_num_cores(data)
                ds_pct = "42." + "{ds_pct:.3}".format(ds_pct=ds_pct).replace("0.", "")
                cmd = ("{samtools} view -O BAM -@ {num_cores} -o {tx_out_file} "
                       "-s {ds_pct} {in_bam}")
                do.run(cmd.format(**locals()), "Downsample BAM file: %s" % os.path.basename(in_bam))
        return out_file

def get_maxcov_downsample_cl(data, in_pipe=None):
    """Retrieve command line for max coverage downsampling, fitting into bamsormadup output.
    """
    max_cov = _get_maxcov_downsample(data) if dd.get_aligner(data) not in ["snap"] else None
    if max_cov:
        if in_pipe == "bamsormadup":
            prefix = "level=0"
        elif in_pipe == "samtools":
            prefix = "-l 0"
        else:
            prefix = ""
        # Swap over to multiple cores until after testing
        #core_arg = "-t %s" % dd.get_num_cores(data)
        core_arg = ""
        return ("%s | variant - -b %s --mark-as-qc-fail --max-coverage %s"
                % (prefix, core_arg, max_cov))
    else:
        if in_pipe == "bamsormadup":
            prefix = "indexfilename={tx_out_file}.bai"
        else:
            prefix = ""
        return prefix

def _get_maxcov_downsample(data):
    """Calculate maximum coverage downsampling for whole genome samples.

    Returns None if we're not doing downsampling.
    """
    from bcbio.bam import ref
    from bcbio.ngsalign import alignprep, bwa
    from bcbio.variation import coverage
    fastq_file = data["files"][0]
    params = alignprep.get_downsample_params(data)
    if params:
        num_reads = alignprep.total_reads_from_grabix(fastq_file)
        if num_reads:
            vrs = dd.get_variant_regions_merged(data)
            total_size = sum([c.size for c in ref.file_contigs(dd.get_ref_file(data), data["config"])])
            if vrs:
                callable_size = pybedtools.BedTool(vrs).total_coverage()
                genome_cov_pct = callable_size / float(total_size)
            else:
                callable_size = total_size
                genome_cov_pct = 1.0
            if (genome_cov_pct > coverage.GENOME_COV_THRESH
                  and dd.get_coverage_interval(data) in ["genome", None, False]):
                total_counts, total_sizes = 0, 0
                for count, size in bwa.fastq_size_output(fastq_file, 5000):
                    total_counts += int(count)
                    total_sizes += (int(size) * int(count))
                read_size = float(total_sizes) / float(total_counts)
                avg_cov = float(num_reads * read_size) / callable_size
                if avg_cov >= params["min_coverage_for_downsampling"]:
                    return int(avg_cov * params["maxcov_downsample_multiplier"])
    return None


def check_header(in_bam, rgnames, ref_file, config):
    """Ensure passed in BAM header matches reference file and read groups names.
    """
    _check_bam_contigs(in_bam, ref_file, config)
    _check_sample(in_bam, rgnames)

def _check_sample(in_bam, rgnames):
    """Ensure input sample name matches expected run group names.
    """
    with pysam.Samfile(in_bam, "rb") as bamfile:
        rg = bamfile.header.get("RG", [{}])
    msgs = []
    warnings = []
    if len(rg) > 1:
        warnings.append("Multiple read groups found in input BAM. Expect single RG per BAM.")
    if len(rg) == 0:
        msgs.append("No read groups found in input BAM. Expect single RG per BAM.")
    if len(rg) > 0 and any(x.get("SM") != rgnames["sample"] for x in rg):
        msgs.append("Read group sample name (SM) does not match configuration `description`: %s vs %s"
                    % (rg[0].get("SM"), rgnames["sample"]))
    if len(msgs) > 0:
        raise ValueError("Problems with pre-aligned input BAM file: %s\n" % (in_bam)
                         + "\n".join(msgs) +
                         "\nSetting `bam_clean: fixrg`\n"
                         "in the configuration can often fix this issue.")
    if warnings:
        print("*** Potential problems in input BAM compared to reference:\n%s\n" %
              "\n".join(warnings))

def _check_bam_contigs(in_bam, ref_file, config):
    """Ensure a pre-aligned BAM file matches the expected reference genome.
    """
    # GATK allows chromosome M to be in multiple locations, skip checking it
    allowed_outoforder = ["chrM", "MT"]
    ref_contigs = [c.name for c in ref.file_contigs(ref_file, config)]
    with pysam.Samfile(in_bam, "rb") as bamfile:
        bam_contigs = [c["SN"] for c in bamfile.header["SQ"]]
    extra_bcs = [x for x in bam_contigs if x not in ref_contigs]
    extra_rcs = [x for x in ref_contigs if x not in bam_contigs]
    problems = []
    warnings = []
    for bc, rc in zip_longest([x for x in bam_contigs if (x not in extra_bcs and
                                                                     x not in allowed_outoforder)],
                                         [x for x in ref_contigs if (x not in extra_rcs and
                                                                     x not in allowed_outoforder)]):
        if bc != rc:
            if bc and rc:
                problems.append("Reference mismatch. BAM: %s Reference: %s" % (bc, rc))
            elif bc:
                warnings.append("Extra BAM chromosomes: %s" % bc)
            elif rc:
                warnings.append("Extra reference chromosomes: %s" % rc)
    for bc in extra_bcs:
        warnings.append("Extra BAM chromosomes: %s" % bc)
    for rc in extra_rcs:
        warnings.append("Extra reference chromosomes: %s" % rc)
    if problems:
        raise ValueError("Unexpected order, name or contig mismatches between input BAM and reference file:\n%s\n"
                         "Setting `bam_clean: remove_extracontigs` in the configuration can often fix this issue."
                         % "\n".join(problems))
    if warnings:
        print("*** Potential problems in input BAM compared to reference:\n%s\n" %
              "\n".join(warnings))

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

def sam_to_bam(in_sam, config):
    if is_bam(in_sam):
        return in_sam

    assert is_sam(in_sam), "%s is not a SAM file" % in_sam
    out_file = os.path.splitext(in_sam)[0] + ".bam"
    if utils.file_exists(out_file):
        return out_file

    samtools = config_utils.get_program("samtools", config)
    num_cores = config["algorithm"].get("num_cores", 1)
    with file_transaction(config, out_file) as tx_out_file:
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
    with file_transaction(config, out_file) as tx_out_file:
        cmd = "{samtools} view -@ {num_cores} -h {in_file} -o {tx_out_file}"
        do.run(cmd.format(**locals()),
               ("Convert BAM to SAM (%s cores): %s to %s"
                % (str(num_cores), in_file, out_file)))
    return out_file


def reheader(header, bam_file, config):
    samtools = config_utils.get_program("samtools", config)
    base, ext = os.path.splitext(bam_file)
    out_file = base + ".reheadered" + ext
    cmd = "{samtools} reheader {header} {bam_file} > {out_file}"
    do.run(cmd.format(**locals()), "Reheadering %s." % bam_file)
    return out_file


def merge(bamfiles, out_bam, config):
    assert all(map(is_bam, bamfiles)), ("Not all of the files to merge are not BAM "
                                        "files: %s " % (bamfiles))
    assert all(map(utils.file_exists, bamfiles)), ("Not all of the files to merge "
                                                   "exist: %s" % (bamfiles))
    if len(bamfiles) == 1:
        return bamfiles[0]
    if os.path.exists(out_bam):
        return out_bam
    samtools = config_utils.get_program("samtools", config)
    bamtools = config_utils.get_program("bamtools", config)
    num_cores = config["algorithm"].get("num_cores", 1)
    with file_transaction(config, out_bam) as tx_out_bam:
        cmd = "{samtools} merge -@ {num_cores} {tx_out_bam} " + " ".join(bamfiles)
        do.run(cmd.format(**locals()), "Merge %s into %s." % (bamfiles, out_bam))
    index(out_bam, config)
    return out_bam


def sort(in_bam, config, order="coordinate", out_dir=None):
    """Sort a BAM file, skipping if already present.
    """
    assert is_bam(in_bam), "%s in not a BAM file" % in_bam
    if bam_already_sorted(in_bam, config, order):
        return in_bam

    sort_stem = _get_sort_stem(in_bam, order, out_dir)
    sort_file = sort_stem + ".bam"
    if not utils.file_exists(sort_file):
        samtools = config_utils.get_program("samtools", config)
        cores = config["algorithm"].get("num_cores", 1)
        with file_transaction(config, sort_file) as tx_sort_file:
            tx_sort_stem = os.path.splitext(tx_sort_file)[0]
            tx_dir = utils.safe_makedir(os.path.dirname(tx_sort_file))
            order_flag = "-n" if order == "queryname" else ""
            resources = config_utils.get_resources("samtools", config)
            # Slightly decrease memory and allow more accurate representation
            # in Mb to ensure fits within systems like SLURM
            mem = config_utils.adjust_memory(resources.get("memory", "2G"),
                                             1.25, "decrease", out_modifier="M").upper()
            cmd = ("{samtools} sort -@ {cores} -m {mem} -O BAM {order_flag} "
                   "-T {tx_sort_stem}-sort -o {tx_sort_file} {in_bam}")
            do.run(cmd.format(**locals()), "Sort BAM file %s: %s to %s" %
                   (order, os.path.basename(in_bam), os.path.basename(sort_file)))
    return sort_file


def bam_already_sorted(in_bam, config, order):
    return order == _get_sort_order(in_bam, config)


def _get_sort_order(in_bam, config):
    for line in pysam.view("-H", in_bam).split("\r\n"):
        if line.startswith("@HD"):
            for keyval in line.split()[1:]:
                key, val = keyval.split(":")
                if key == "SO":
                    return val

def _get_sort_stem(in_bam, order, out_dir):
    SUFFIXES = {"coordinate": ".sorted", "queryname": ".nsorted"}
    sort_base = os.path.splitext(in_bam)[0]
    if out_dir:
        sort_base = os.path.join(out_dir, os.path.basename(sort_base))
    for suffix in SUFFIXES:
        sort_base = sort_base.split(suffix)[0]
    return sort_base + SUFFIXES[order]

def aligner_from_header(in_bam):
    """Identify aligner from the BAM header; handling pre-aligned inputs.
    """
    from bcbio.pipeline.alignment import TOOLS
    with pysam.Samfile(in_bam, "rb") as bamfile:
        for pg in bamfile.header.get("PG", []):
            for ka in TOOLS.keys():
                if pg.get("PN", "").lower().find(ka) >= 0:
                    return ka

def sample_name(in_bam):
    """Get sample name from BAM file.
    """
    with pysam.AlignmentFile(in_bam, "rb", check_sq=False) as in_pysam:
        try:
            if "RG" in in_pysam.header:
                return in_pysam.header["RG"][0]["SM"]
        except ValueError:
            return None

def estimate_read_length(bam_file, nreads=1000):
    """
    estimate median read length of a SAM/BAM file
    """
    with open_samfile(bam_file) as bam_handle:
        reads = tz.itertoolz.take(nreads, bam_handle)
        lengths = [len(x.seq) for x in reads]
    return int(numpy.median(lengths))

def estimate_fragment_size(bam_file, nreads=5000):
    """
    estimate median fragment size of a SAM/BAM file
    """
    with open_samfile(bam_file) as bam_handle:
        reads = tz.itertoolz.take(nreads, bam_handle)
        # it would be good to skip spliced paired reads.
        lengths = [x.template_length for x in reads if x.template_length > 0]
    if not lengths:
        return 0
    return int(numpy.median(lengths))

def filter_primary(bam_file, data):
    """Filter reads to primary only BAM.

    Removes:
      - not primary alignment (0x100) 256
      - supplementary alignment (0x800) 2048
    """
    stem, ext = os.path.splitext(bam_file)
    out_file = stem + ".primary" + ext
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cores = dd.get_num_cores(data)
            cmd = ("samtools view -@ {cores} -F 2304 -b {bam_file} > {tx_out_file}")
            do.run(cmd.format(**locals()), ("Filtering primary alignments in %s." %
                                            os.path.basename(bam_file)))
    return out_file

def estimate_max_mapq(in_bam, nreads=1e6):
    """Guess maximum MAPQ in a BAM file of reads with alignments
    """
    with pysam.Samfile(in_bam, "rb") as work_bam:
        reads = tz.take(int(nreads), work_bam)
        return max([x.mapq for x in reads if not x.is_unmapped])

def convert_cufflinks_mapq(in_bam, out_bam=None):
    """Cufflinks expects the not-valid 255 MAPQ for uniquely mapped reads.
    This detects the maximum mapping quality in a BAM file and sets
    reads with that quality to be 255
    """
    CUFFLINKSMAPQ = 255
    if not out_bam:
        out_bam = os.path.splitext(in_bam)[0] + "-cufflinks.bam"
    if utils.file_exists(out_bam):
        return out_bam
    maxmapq = estimate_max_mapq(in_bam)
    if maxmapq == CUFFLINKSMAPQ:
        return in_bam
    logger.info("Converting MAPQ scores in %s to be Cufflinks compatible." % in_bam)
    with pysam.Samfile(in_bam, "rb") as in_bam_fh:
        with pysam.Samfile(out_bam, "wb", template=in_bam_fh) as out_bam_fh:
            for read in in_bam_fh:
                if read.mapq == maxmapq and not read.is_unmapped:
                    read.mapq = CUFFLINKSMAPQ
                out_bam_fh.write(read)
    return out_bam

def convert_invalid_mapq(in_bam, out_bam=None):
    """Some aligners output 255 to denote a uniquely mapped read which is
    an invalid MAPQ value according to the SAM spec. This detects
    that and changes it to be 60.
    """
    INVALIDMAPQ = 255
    VALIDMAPQ = 60
    if not out_bam:
        out_bam = os.path.splitext(in_bam)[0] + "-MAPQfixed.bam"
    if utils.file_exists(out_bam):
        return out_bam
    maxmapq = estimate_max_mapq(in_bam)
    if maxmapq != INVALIDMAPQ:
        return in_bam
    logger.info("Converting 255 MAPQ scores in %s to 60." % in_bam)
    with pysam.Samfile(in_bam, "rb") as in_bam_fh:
        with pysam.Samfile(out_bam, "wb", template=in_bam_fh) as out_bam_fh:
            for read in in_bam_fh:
                if read.mapq == INVALIDMAPQ and not read.is_unmapped:
                    read.mapq = VALIDMAPQ
                out_bam_fh.write(read)
    return out_bam
