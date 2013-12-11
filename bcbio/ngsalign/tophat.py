"""Next-gen alignments with TopHat a spliced read mapper for RNA-seq experiments.

http://tophat.cbcb.umd.edu
"""
import sh
import os
import shutil
from contextlib import closing
import glob

from py_descriptive_statistics import Enum as Stats
import pysam

from bcbio.pipeline import config_utils
from bcbio.ngsalign import bowtie, bowtie2
from bcbio.utils import safe_makedir, file_exists, get_in
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger

from bcbio.provenance import do
from bcbio import broad
from bcbio.broad.metrics import PicardMetricsParser

_out_fnames = ["accepted_hits.sam", "junctions.bed",
               "insertions.bed", "deletions.bed"]


def _set_quality_flag(options, config):
    qual_format = config["algorithm"].get("quality_format", None)
    if qual_format.lower() == "illumina":
        options["solexa1.3-quals"] = True
    elif qual_format.lower() == "solexa":
        options["solexa-quals"] = True
    return options

def _set_transcriptome_option(options, data, ref_file):
    # prefer transcriptome-index vs a GTF file if available
    transcriptome_index = get_in(data, ("genome_resources", "rnaseq",
                                        "transcriptome_index", "tophat"))
    if transcriptome_index and file_exists(transcriptome_index):
        options["transcriptome-index"] = os.path.splitext(transcriptome_index)[0]
        return options

    gtf_file = data["genome_resources"]["rnaseq"].get("transcripts")
    if gtf_file:
        options["GTF"] = gtf_file
        return options

    return options

def _set_cores(options, config):
    num_cores = config["algorithm"].get("num_cores", 0)
    if num_cores > 1 and "num-threads" not in options:
        options["num-threads"] = num_cores
    return options

def _set_rg_options(options, names):
    if not names:
        return options
    options["rg-id"] = names["rg"]
    options["rg-sample"] = names["sample"]
    options["rg-library"] = names["pl"]
    options["rg-platform-unit"] = names["pu"]
    return options

def _set_stranded_flag(options, config):
    strand_flag = {"unstranded": "fr-unstranded",
                   "firststrand": "fr-firststrand",
                   "secondstrand": "fr-secondstrand"}
    stranded = get_in(config, ("algorithm", "strandedness"), "unstranded").lower()
    assert stranded in strand_flag, ("%s is not a valid strandedness value. "
                                     "Valid values are 'firststrand', "
                                     "'secondstrand' and 'unstranded" % (stranded))
    flag = strand_flag[stranded]
    options["library-type"] = flag
    return options

def tophat_align(fastq_file, pair_file, ref_file, out_base, align_dir, data,
                 names=None):
    """
    run alignment using Tophat v2
    """
    config = data["config"]
    options = get_in(config, ("resources", "tophat", "options"), {})
    options = _set_quality_flag(options, config)
    options = _set_transcriptome_option(options, data, ref_file)
    options = _set_cores(options, config)
    options = _set_rg_options(options, names)
    options = _set_stranded_flag(options, config)

    # select the correct bowtie option to use; tophat2 is ignoring this option
    if _tophat_major_version(config) == 1:
        raise NotImplementedError("Tophat versions < 2.0 are not supported, please "
                                  "download the newest version of Tophat here: "
                                  "http://tophat.cbcb.umd.edu")
    if _ref_version(ref_file) == 1:
        options["bowtie1"] = True

    out_dir = os.path.join(align_dir, "%s_tophat" % out_base)
    out_file = os.path.join(out_dir, _out_fnames[0])
    files = [ref_file, fastq_file]
    if not file_exists(out_file):
        with file_transaction(out_dir) as tx_out_dir:
            _check_bowtie(ref_file, config)
            safe_makedir(tx_out_dir)
            if pair_file and not options.get("mate-inner-dist", None):
                d, d_stdev = _estimate_paired_innerdist(fastq_file, pair_file,
                                                        ref_file, out_base,
                                                        tx_out_dir, data)
                options["mate-inner-dist"] = d
                options["mate-std-dev"] = d_stdev
                files.append(pair_file)
            options["output-dir"] = tx_out_dir
            options["no-convert-bam"] = True
            options["no-coverage-search"] = True
            tophat_runner = sh.Command(config_utils.get_program("tophat",
                                                                config))
            ready_options = {}
            for k, v in options.iteritems():
                ready_options[k.replace("-", "_")] = v
            # tophat requires options before arguments,
            # otherwise it silently ignores them
            tophat_ready = tophat_runner.bake(**ready_options)
            cmd = str(tophat_ready.bake(*files))
            do.run(cmd, "Running Tophat on %s and %s." % (fastq_file, pair_file), None)
        _fix_empty_readnames(out_file)
    if pair_file and _has_alignments(out_file):
        final_out = _fix_mates(out_file, os.path.join(out_dir, "%s-align.bam" % out_base),
                               ref_file, config)
    else:
        final_out = os.path.join(out_dir, "%s.sam" % out_base)
        if not file_exists(final_out):
            os.symlink(os.path.basename(out_file), final_out)
    return final_out

def _has_alignments(sam_file):
    print sam_file
    with open(sam_file) as in_handle:
        for line in in_handle:
            if line.startswith("File removed to save disk space"):
                return False
            elif not line.startswith("@"):
                return True
    return False

def _fix_empty_readnames(orig_file):
    """ Fix SAMfile reads with empty read names

    Tophat 2.0.9 sometimes outputs empty read names, making the
    FLAG field be the read name. This throws those reads away.
    """
    with file_transaction(orig_file) as tx_out_file:
        logger.info("Removing reads with empty read names from Tophat output.")
        with open(orig_file) as orig, open(tx_out_file, "w") as out:
            for line in orig:
                if line.split()[0].isdigit():
                    continue
                out.write(line)
    return orig_file


def _fix_mates(orig_file, out_file, ref_file, config):
    """Fix problematic unmapped mate pairs in TopHat output.

    TopHat 2.0.9 appears to have issues with secondary reads:
    https://groups.google.com/forum/#!topic/tuxedo-tools-users/puLfDNbN9bo
    This cleans the input file to only keep properly mapped pairs,
    providing a general fix that will handle correctly mapped secondary
    reads as well.
    """
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            samtools = config_utils.get_program("samtools", config)
            cmd = "{samtools} view -bt {ref_file}.fai -F 8 {orig_file} > {tx_out_file}"
            do.run(cmd.format(**locals()), "Fix mate pairs in TopHat output", {})
    return out_file

def align(fastq_file, pair_file, ref_file, out_base, align_dir, data,
          names=None):
    out_files = tophat_align(fastq_file, pair_file, ref_file, out_base,
                             align_dir, data, names)

    return out_files


def _estimate_paired_innerdist(fastq_file, pair_file, ref_file, out_base,
                               out_dir, data):
    """Use Bowtie to estimate the inner distance of paired reads.
    """
    # skip initial reads for large file, but not for smaller
    # mean, stdev = _bowtie_for_innerdist("1000000", fastq_file, pair_file, ref_file,
    #                              out_base, out_dir, config)
    # if it is a small file, use the old method
    mean, stdev = _small_file_innerdist("100000", fastq_file, pair_file, ref_file,
                                        out_base, out_dir, data, True)
    if not mean or not stdev:
        mean, stdev = _small_file_innerdist("1", fastq_file, pair_file, ref_file,
                                            out_base, out_dir, data, True)
    # No reads aligning so no data to process, set some default values
    if not mean or not stdev:
        mean, stdev = 200, 50

    return mean, stdev


def _bowtie_for_innerdist(start, fastq_file, pair_file, ref_file, out_base,
                          out_dir, config, remove_workdir=False):
    work_dir = os.path.join(out_dir, "innerdist_estimate")
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    safe_makedir(work_dir)
    extra_args = ["-s", str(start), "-u", "250000"]
    bowtie_runner = _select_bowtie_version(config)
    out_sam = bowtie_runner.align(fastq_file, pair_file, ref_file, out_base,
                                  work_dir, config, extra_args)
    runner = broad.runner_from_config(config)
    metrics_file = runner.run_fn("picard_insert_metrics", out_sam)
    if not file_exists(metrics_file):
        return None, None
    parser = PicardMetricsParser()
    with open(metrics_file) as metrics_handle:
        insert_metrics = parser._parse_insert_metrics(metrics_handle)

    avg_read_length = _calculate_average_read_length(out_sam)
    mean_insert = int(float(insert_metrics["MEAN_INSERT_SIZE"])) - int(2 * avg_read_length)
    std_deviation = int(float(insert_metrics["STANDARD_DEVIATION"]))
    return mean_insert, std_deviation

def _small_file_innerdist(start, fastq_file, pair_file, ref_file, out_base,
                          out_dir, data, remove_workdir=False):
    work_dir = os.path.join(out_dir, "innerdist_estimate")
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    safe_makedir(work_dir)
    extra_args = ["-s", str(start), "-u", "250000"]
    bowtie_runner = _select_bowtie_version(data["config"])
    out_sam = bowtie_runner.align(fastq_file, pair_file, ref_file, out_base,
                                  work_dir, data, extra_args)
    dists = []
    with closing(pysam.Samfile(out_sam)) as work_sam:
        for read in work_sam:
            if read.is_proper_pair and read.is_read1:
                dists.append(abs(read.isize) - 2 * read.rlen)
    if dists:
        dist_stats = Stats(dists)
        return int(round(dist_stats.mean())), int(round(dist_stats.standard_deviation()))
    else:
        return None, None

def _calculate_average_read_length(sam_file):
    with closing(pysam.Samfile(sam_file)) as work_sam:
        count = 0
        read_lengths = []
        for read in work_sam:
            count = count + 1
            read_lengths.append(read.rlen)
    avg_read_length = int(float(sum(read_lengths)) / float(count))
    return avg_read_length


def _check_bowtie(ref_file, config):
    if not _bowtie_ref_match(ref_file, config):
        logger.error("Bowtie version %d was detected but the reference "
                     "file %s is built for version %d. Download version "
                     "%d or build it with bowtie-build."
                     % (_bowtie_major_version(config), ref_file,
                        _ref_version(ref_file),
                        _bowtie_major_version(config)))
        exit(1)

def _bowtie_major_version(config):
    bowtie_runner = sh.Command(config_utils.get_program("bowtie", config,
                                                        default="bowtie2"))
    """
    bowtie --version returns strings like this:
    bowtie version 0.12.7
    32-bit
    Built on Franklin.local
    Tue Sep  7 14:25:02 PDT 2010
    """
    version_line = str(bowtie_runner(version=True)).split("\n")[0]
    version_string = version_line.strip().split()[2]
    major_version = int(version_string.split(".")[0])
    # bowtie version 1 has a leading character of 0
    if major_version == 0:
        major_version += 1
    return major_version


def _tophat_major_version(config):
    tophat_runner = sh.Command(config_utils.get_program("tophat", config,
                                                        default="tophat"))

    # tophat --version returns strings like this: Tophat v2.0.4
    version_string = str(tophat_runner(version=True)).strip().split()[1]
    major_version = int(version_string.split(".")[0][1:])
    return major_version


def _bowtie_ref_match(ref_file, config):
    return _ref_version(ref_file) == _bowtie_major_version(config)


def _select_bowtie_version(config):
    if _bowtie_major_version(config) == 1:
        return bowtie
    else:
        return bowtie2


def _ref_version(ref_file):
    for ext in [os.path.splitext(x)[1] for x in glob.glob(ref_file + "*")]:
        if ext == ".ebwt":
            return 1
        elif ext == ".bt2":
            return 2
    raise ValueError("Cannot detect which reference version %s is. "
                     "Should end in either .ebwt (bowtie) or .bt2 "
                     "(bowtie2)." % (ref_file))


def job_requirements(cores, memory):
    MIN_TOPHAT_MEMORY = 8.0
    if not memory or cores * memory < MIN_TOPHAT_MEMORY:
        memory = MIN_TOPHAT_MEMORY / cores
    return cores, memory

align.job_requirements = job_requirements
