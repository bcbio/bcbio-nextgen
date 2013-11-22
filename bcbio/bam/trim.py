"""Provide trimming of input reads from Fastq or BAM files.
"""
import os
import contextlib
from bcbio.utils import (file_exists, save_diskspace, safe_makedir,
                         replace_suffix, append_stem, is_pair,
                         replace_directory, map_wrap)
from bcbio.log import logger
from bcbio.bam import fastq
from bcbio.provenance import do
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from itertools import izip, repeat
from bcbio.distributed.transaction import file_transaction


SUPPORTED_ADAPTERS = {
    "illumina": ["AACACTCTTTCCCT", "AGATCGGAAGAGCG"],
    "truseq": ["AGATCGGAAGAG"],
    "polya": ["AAAAAAAAAAAAA"],
    "nextera": ["AATGATACGGCGA", "CAAGCAGAAGACG"]}

def brun_trim_fastq(fastq_files, dirs, config):
    """Trim FASTQ files, removing low quality B-runs.

    This removes stretches of low quality sequence from read ends. Illumina
    quality assessment generates these stretches. Removing them can help reduce
    false positive rates for variant calling.

    http://genomebiology.com/2011/12/11/R112

    Does simple trimming of problem ends and removes read pairs where
    any of the trimmed read sizes falls below the allowable size.
    """
    qual_format = config["algorithm"].get("quality_format", "").lower()
    min_length = int(config["algorithm"].get("min_read_length", 20))
    to_trim = "B" if qual_format == "illumina" else "#"
    with _work_handles(fastq_files, dirs, "-qtrim.txt") as (in_handles, out_handles, out_fnames):
        if len(out_handles) == len(fastq_files):
            for next_reads in _trim_by_read(in_handles, to_trim, min_length):
                for fname, (name, seq, qual) in next_reads.iteritems():
                    out_handles[fname].write("@%s\n%s\n+\n%s\n" % (name, seq, qual))
        out_files = [out_fnames[x] for x in fastq_files]
        for inf, outf in zip(fastq_files, out_files):
            _save_diskspace(inf, outf, config)
        return out_files

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


def _trim_quality(seq, qual, to_trim, min_length):
    """Trim bases of the given quality from 3' read ends.
    """
    removed = 0
    while qual.endswith(to_trim):
        removed += 1
        qual = qual[:-1]
    if len(qual) >= min_length:
        return seq[:len(seq) - removed], qual
    else:
        return None, None

@contextlib.contextmanager
def _work_handles(in_files, dirs, ext):
    """Create working handles for input files and close on completion.
    """
    out_dir = safe_makedir(os.path.join(dirs["work"], "trim"))
    out_handles = {}
    in_handles = {}
    name_map = {}
    for in_file in in_files:
        out_file = os.path.join(out_dir, "{base}{ext}".format(
            base=os.path.splitext(os.path.basename(in_file))[0], ext=ext))
        name_map[in_file] = out_file
        if not file_exists(out_file):
            in_handles[in_file] = open(in_file)
            out_handles[in_file] = open(out_file, "w")
    try:
        yield in_handles, out_handles, name_map
    finally:
        for h in in_handles.values():
            h.close()
        for h in out_handles.values():
            h.close()

def _trim_by_read(in_handles, to_trim, min_length):
    """Lazy generator for trimmed reads for all input files.
    """
    iterators = [(f, FastqGeneralIterator(h)) for f, h in in_handles.iteritems()]
    f1, x1 = iterators[0]
    for name, seq, qual in x1:
        out = {}
        tseq, tqual = _trim_quality(seq, qual, to_trim, min_length)
        if tseq:
            out[f1] = (name, tseq, tqual)
        for f2, x2 in iterators[1:]:
            name, seq, qual = x2.next()
            tseq, tqual = _trim_quality(seq, qual, to_trim, min_length)
            if tseq:
                out[f2] = (name, tseq, tqual)
        if len(out) == len(iterators):
            yield out

def _save_diskspace(in_file, out_file, config):
    """Potentially remove input file to save space if configured and in work directory.
    """
    if (os.path.commonprefix([in_file, out_file]).rstrip("/") ==
        os.path.split(os.path.dirname(out_file))[0]):
        save_diskspace(in_file, "Trimmed to {}".format(out_file), config)



def remove_short_reads(fastq_files, dirs, lane_config):
    """
    remove reads from a single or pair of fastq files which fall below
    a length threshold (30 bases)

    """
    MIN_LENGTH = 20
    supplied_quality_format = _get_quality_format(lane_config)
    if supplied_quality_format == "illumina":
        quality_format = "fastq-illumina"
    else:
        quality_format = "fastq-sanger"

    if is_pair(fastq_files):
        fastq1, fastq2 = fastq_files
        out_files = fastq.filter_reads_by_length(fastq1, fastq2, quality_format,
                                                 MIN_LENGTH)
    else:
        out_files = [fastq.filter_single_reads_by_length(fastq_files[0],
                                                        quality_format,
                                                        MIN_LENGTH)]
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
    if quality_format == "illumina":
        quality_base = "64"
    else:
        quality_base = "33"

    # --times=2 tries twice remove adapters which will allow things like:
    # realsequenceAAAAAAadapter to remove both the poly-A and the adapter
    # this behavior might not be what we want; we could also do two or
    # more passes of cutadapt
    base_cmd = ["cutadapt", "--times=" + "2", "--quality-base=" + quality_base,
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
