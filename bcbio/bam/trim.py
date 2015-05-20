"""Provide trimming of input reads from Fastq or BAM files.
"""
import os
import sys

from bcbio import utils
from bcbio.utils import (file_exists, append_stem, replace_directory)
from bcbio.log import logger
from bcbio.distributed import objectstore
from bcbio.provenance import do
from Bio.Seq import Seq
from bcbio.distributed.transaction import file_transaction

MINIMUM_LENGTH = 25

SUPPORTED_ADAPTERS = {
    "illumina": ["AACACTCTTTCCCT", "AGATCGGAAGAGCG"],
    "truseq": ["AGATCGGAAGAG"],
    "polya": ["AAAAAAAAAAAAA"],
    "nextera": ["AATGATACGGCGA", "CAAGCAGAAGACG"]}

def trim_adapters(fastq_files, out_dir, config):
    """
    for small insert sizes, the read length can be longer than the insert
    resulting in the reverse complement of the 3' adapter being sequenced.
    this takes adapter sequences and trims the only the reverse complement
    of the adapter

    MYSEQUENCEAAAARETPADA -> MYSEQUENCEAAAA (no polyA trim)

    """
    quality_format = _get_quality_format(config)
    to_trim = _get_sequences_to_trim(config, SUPPORTED_ADAPTERS)
    out_files = replace_directory(append_stem(fastq_files, ".trimmed"), out_dir)
    out_files = _cutadapt_trim(fastq_files, quality_format, to_trim, out_files, config)
    return out_files

def _cutadapt_trim(fastq_files, quality_format, adapters, out_files, config):
    """Trimming with cutadapt, using version installed with bcbio-nextgen.

    Uses the system executable to find the version next to our Anaconda Python.
    TODO: Could we use cutadapt as a library to avoid this?
    """
    if all([file_exists(x) for x in out_files]):
        return out_files
    cmd = _cutadapt_trim_cmd(fastq_files, quality_format, adapters, out_files)
    if len(fastq_files) == 1:
        of1 = out_files[0]
        message = "Trimming %s in single end mode with cutadapt." % (fastq_files[0])
        with file_transaction(config, of1) as of1_tx:
            do.run(cmd.format(**locals()), message)
    else:
        with file_transaction(config, out_files) as tx_out_files:
            of1_tx, of2_tx = tx_out_files
            tmp_fq1 = append_stem(of1_tx, ".tmp")
            tmp_fq2 = append_stem(of2_tx, ".tmp")
            singles_file = of1_tx + ".single"
            message = "Trimming %s and %s in paired end mode with cutadapt." % (fastq_files[0],
                                                                                fastq_files[1])
            do.run(cmd.format(**locals()), message)
    return out_files

def _cutadapt_trim_cmd(fastq_files, quality_format, adapters, out_files):
    """Trimming with cutadapt, using version installed with bcbio-nextgen.

    Uses the system executable to find the version next to our Anaconda Python.
    TODO: Could we use cutadapt as a library to avoid this?
    """
    if all([file_exists(x) for x in out_files]):
        return out_files
    if quality_format == "illumina":
        quality_base = "64"
    else:
        quality_base = "33"

    # --times=2 tries twice remove adapters which will allow things like:
    # realsequenceAAAAAAadapter to remove both the poly-A and the adapter
    # this behavior might not be what we want; we could also do two or
    # more passes of cutadapt
    cutadapt = os.path.join(os.path.dirname(sys.executable), "cutadapt")
    adapter_cmd = " ".join(map(lambda x: "--adapter=" + x, adapters))
    base_cmd = ("{cutadapt} --times=2 --quality-base={quality_base} "
                "--quality-cutoff=5 --format=fastq "
                "{adapter_cmd} ").format(**locals())
    if len(fastq_files) == 1:
        return _cutadapt_se_cmd(fastq_files, out_files, base_cmd)
    else:
        return _cutadapt_pe_nosickle(fastq_files, out_files, quality_format, base_cmd)

def _cutadapt_se_cmd(fastq_files, out_files, base_cmd):
    """
    this has to use the -o option, not redirect to stdout in order for gzipping to be
    honored
    """
    min_length = MINIMUM_LENGTH
    cmd = base_cmd + " --minimum-length={min_length} ".format(**locals())
    fq1 = objectstore.cl_input(fastq_files[0])
    of1 = out_files[0]
    cmd += " -o {of1} " + str(fq1)
    return cmd

def _cutadapt_pe_nosickle(fastq_files, out_files, quality_format, base_cmd):
    """
    sickle has an issue with 0 length reads, here is the open issue for it:
    https://github.com/najoshi/sickle/issues/32
    until that is resolved, this is a workaround which avoids using sickle
    """
    fq1, fq2 = [objectstore.cl_input(x) for x in fastq_files]
    of1, of2 = out_files
    base_cmd += " --minimum-length={min_length} ".format(min_length=MINIMUM_LENGTH)
    first_cmd = base_cmd + " -o {tmp_fq1} -p {tmp_fq2} " + fq1 + " " + fq2
    second_cmd = base_cmd + " -o {of2_tx} -p {of1_tx} {tmp_fq2} {tmp_fq1}"
    return first_cmd + ";" + second_cmd + "; rm {tmp_fq1} {tmp_fq2} "

def _get_sequences_to_trim(config, builtin):
    builtin_adapters = _get_builtin_adapters(config, builtin)
    polya = builtin_adapters.get("polya", [None])[0]
    # allow for trimming of custom sequences for advanced users
    custom_trim = config["algorithm"].get("custom_trim", [])
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

def _get_quality_format(config):
    SUPPORTED_FORMATS = ["illumina", "standard"]
    quality_format = config["algorithm"].get("quality_format", "standard").lower()
    if quality_format not in SUPPORTED_FORMATS:
        logger.error("quality_format is set to an unsupported format. "
                     "Supported formats are %s."
                     % (", ".join(SUPPORTED_FORMATS)))
        exit(1)
    return quality_format

def _get_builtin_adapters(config, builtin):
    chemistries = config["algorithm"].get("adapters", [])
    adapters = {chemistry: builtin[chemistry] for
                chemistry in chemistries if chemistry in builtin}
    return adapters
