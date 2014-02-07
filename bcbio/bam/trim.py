"""Provide trimming of input reads from Fastq or BAM files.
"""
import os
import sys

from bcbio.utils import (file_exists, safe_makedir,
                         replace_suffix, append_stem, is_pair,
                         replace_directory, map_wrap)
from bcbio.log import logger
from bcbio.bam import fastq
from bcbio.provenance import do
from Bio.Seq import Seq
from itertools import izip, repeat
from bcbio.distributed.transaction import file_transaction


SUPPORTED_ADAPTERS = {
    "illumina": ["AACACTCTTTCCCT", "AGATCGGAAGAGCG"],
    "truseq": ["AGATCGGAAGAG"],
    "polya": ["AAAAAAAAAAAAA"],
    "nextera": ["AATGATACGGCGA", "CAAGCAGAAGACG"]}

def trim_read_through(fastq_files, dirs, lane_config):
    """
    for small insert sizes, the read length can be longer than the insert
    resulting in the reverse complement of the 3' adapter being sequenced.
    this takes adapter sequences and trims the only the reverse complement
    of the adapter

    MYSEQUENCEAAAARETPADA -> MYSEQUENCEAAAA (no polyA trim)

    """
    quality_format = _get_quality_format(lane_config)
    to_trim = _get_sequences_to_trim(lane_config)
    out_files = _get_read_through_trimmed_outfiles(fastq_files, dirs)
    fixed_files = append_stem(out_files, ".fixed")
    if all(map(file_exists, fixed_files)):
        return fixed_files
    logger.info("Trimming %s from the 3' end of reads in %s using "
                "cutadapt." % (", ".join(to_trim),
                               ", ".join(fastq_files)))
    cores = lane_config["algorithm"].get("num_cores", 1)
    out_files = _cutadapt_trim(fastq_files, quality_format,
                               to_trim, out_files, cores)

    fixed_files = remove_short_reads(out_files, dirs, lane_config)
    return fixed_files

def remove_short_reads(fastq_files, dirs, lane_config):
    """
    remove reads from a single or pair of fastq files which fall below
    a length threshold (30 bases)

    """
    min_length = int(lane_config["algorithm"].get("min_read_length", 20))
    supplied_quality_format = _get_quality_format(lane_config)
    if supplied_quality_format == "illumina":
        quality_format = "fastq-illumina"
    else:
        quality_format = "fastq-sanger"

    if is_pair(fastq_files):
        fastq1, fastq2 = fastq_files
        out_files = fastq.filter_reads_by_length(fastq1, fastq2, quality_format, min_length)
    else:
        out_files = [fastq.filter_single_reads_by_length(fastq_files[0],
                                                         quality_format, min_length)]
    map(os.remove, fastq_files)
    return out_files

def _get_read_through_trimmed_outfiles(fastq_files, dirs):
    out_dir = os.path.join(dirs["work"], "trim")
    safe_makedir(out_dir)
    out_files = replace_directory(append_stem(fastq_files, "_trimmed"),
                                  out_dir)
    return out_files

def _get_sequences_to_trim(lane_config):
    builtin_adapters = _get_builtin_adapters(lane_config)
    polya = builtin_adapters.get("polya", [None])[0]
    # allow for trimming of custom sequences for advanced users
    custom_trim = lane_config["algorithm"].get("custom_trim", [])
    builtin_adapters = {k: v for k, v in builtin_adapters.items() if
                        k != "polya"}
    trim_sequences = custom_trim
    # for unstranded RNA-seq, libraries, both polyA and polyT can appear
    # at the 3' end as well
    if polya:
        trim_sequences += [polya, str(Seq(polya).reverse_complement())]

    # also trim the reverse complement of the adapters
    for _, v in builtin_adapters.items():
        trim_sequences += [str(Seq(sequence)) for sequence in v]
        trim_sequences += [str(Seq(sequence).reverse_complement()) for
                           sequence in v]
    return trim_sequences


def _cutadapt_trim(fastq_files, quality_format, adapters, out_files, cores):
    """Trimming with cutadapt, using version installed with bcbio-nextgen.

    Uses the system executable to find the version next to our Anaconda Python.
    TODO: Could we use cutadapt as a library to avoid this?
    """
    if quality_format == "illumina":
        quality_base = "64"
    else:
        quality_base = "33"

    # --times=2 tries twice remove adapters which will allow things like:
    # realsequenceAAAAAAadapter to remove both the poly-A and the adapter
    # this behavior might not be what we want; we could also do two or
    # more passes of cutadapt
    cutadapt = os.path.join(os.path.dirname(sys.executable), "cutadapt")
    base_cmd = [cutadapt, "--times=" + "2", "--quality-base=" + quality_base,
                "--quality-cutoff=20", "--format=fastq", "--minimum-length=0"]
    adapter_cmd = map(lambda x: "--adapter=" + x, adapters)
    base_cmd.extend(adapter_cmd)
    if all(map(file_exists, out_files)):
        return out_files
    with file_transaction(out_files) as tmp_out_files:
        if isinstance(tmp_out_files, basestring):
            tmp_out_files = [tmp_out_files]
        map(_run_cutadapt_on_single_file, izip(repeat(base_cmd), fastq_files,
                                               tmp_out_files))
    return out_files

@map_wrap
def _run_cutadapt_on_single_file(base_cmd, fastq_file, out_file):
    stat_file = replace_suffix(out_file, ".trim_stats.txt")
    with open(stat_file, "w") as stat_handle:
        cmd = list(base_cmd)
        cmd.extend(["--output=" + out_file, fastq_file])
        do.run(cmd, "Running cutadapt on %s." % (fastq_file), None)


def _get_quality_format(lane_config):
    SUPPORTED_FORMATS = ["illumina", "standard"]
    quality_format = lane_config["algorithm"].get("quality_format",
                                                  "standard").lower()
    if quality_format not in SUPPORTED_FORMATS:
        logger.error("quality_format is set to an unsupported format. "
                     "Supported formats are %s."
                     % (", ".join(SUPPORTED_FORMATS)))
        exit(1)
    return quality_format

def _get_builtin_adapters(lane_config):
    chemistries = lane_config["algorithm"].get("adapters", [])
    adapters = {chemistry: SUPPORTED_ADAPTERS[chemistry] for
                chemistry in chemistries if chemistry in SUPPORTED_ADAPTERS}
    return adapters
