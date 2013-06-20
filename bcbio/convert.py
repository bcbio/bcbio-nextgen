"""
file format conversion utility functions
"""
import pysam
from bcbio.utils import (file_exists, replace_suffix, tmpfile)
from bcbio.distributed.transaction import file_transaction
import os
import subprocess
from bcbio.log import logger
from bcbio.pipeline.config_utils import get_transcript_gtf

def expects(fmt):
    """
    decorates a function with a file format checker. for example if your function
    foo only excepts bam files:

    @expects("bam")
    def foo(in_file):
        return "bam"

    foo("in.bam") -> "bam"
    foo("in.sam") -> ValueError exception thrown

    requires two things:
    1) in_file is the first argument or your function takes a keyword argument in_file
    2) you have implemented an is_fmt predicate (for example is_bam) that returns True
    if in_file is the right format
    """
    def decorator(fn):
        def wrapped(*args, **kwargs):
            in_file = kwargs.get("in_file", None)
            if not in_file:
                in_file = args[0]
            id_fun_string = "is_" + fmt
            id_fun = globals().get(id_fun_string)
            if not id_fun:
                _format_detector_not_found(id_fun_string)
            if not id_fun(in_file):
                _wrong_input_file_format(fn, fmt, in_file)
            return fn(*args, **kwargs)
        return wrapped
    return decorator

def _format_detector_not_found(id_fun_string):
    error_string = "The file format detector {0} does not exist."
    raise ValueError(error_string.format(id_fun_string))

def _wrong_input_file_format(fn, fmt, in_file):
    error_string = "{0} expects {1} but {2} is not of that format."
    raise ValueError(error_string.format(fn, fmt, in_file))

def bam2sam(in_file, samtools="samtools"):
    """
    converts a bam file to a sam file
    bam2sam("file.bam") -> "file.sam"
    """
    assert(is_bam(in_file)), "bam2sam requires a BAM file, got %s" % in_file
    out_file = replace_suffix(in_file, ".sam")
    if file_exists(out_file):
        return out_file
    with file_transaction(out_file) as tmp_out_file:
        #pysam.view("-h", "-o" + tmp_out_file, in_file)
        cmd = "{samtools} view -h -o {tmp_out_file} {in_file}".format(**locals())
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError:
            cmd_string = subprocess.list2cmdline(cmd)
            logger.error("bam2sam returned an error. The command "
                         "used to run bam2sam was: %s." % (cmd_string))
    return out_file

def bamindex(in_file, samtools="samtools"):
    """
    index a bam file
    avoids use of pysam.index which is not working for indexing as of 0.7.4
    with ipython
    """
    assert(is_bam(in_file)), "bamindex requires a BAM file, got %s" % in_file
    out_file = replace_suffix(in_file, ".bai")
    if file_exists(out_file):
        return out_file
    cmd = ["samtools", "index", in_file]
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError:
        cmd_string = subprocess.list2cmdline(cmd)
        logger.error("bamindex returned an error. The command "
                     "used to run bamindex was: %s." % (cmd_string))
    return out_file

def bam2sizes(in_file):
    """
    converts a bam file to a chromosome sizes file
    chromosome sizes has the format:
    CHROM	SIZE

    For example:
    chr1	249250621
    chr2	243199373

    bam2sizes("in.bam") -> "in.sizes"
    """
    assert(is_bam(in_file)), "bam2sizes requires a BAM file, got %s" % in_file
    base, _ = os.path.splitext(in_file)
    out_file = base + ".sizes"
    if file_exists(out_file):
        return out_file

    header = pysam.Samfile(in_file, 'rb').header
    with file_transaction(out_file) as tmp_out_file:
        with open(tmp_out_file, 'w') as out_handle:
            for line in header['SQ']:
                out_handle.write("\t".join([line['SN'], str(line['LN'])]) + "\n")
    return out_file

def bam2freec_len(in_file):
    """
    converts a bam file to a chromosome length file for the use with
    Control-FREEC:
    http://bioinfo-out.curie.fr/projects/freec/
    the length file has the format:
    INDEX	CHROM	SIZE

    For example:
    1	chr1	249250621
    2	chr2	243199373

    bam2freec_len("in.bam") -> "in.len"
    """
    assert(is_bam(in_file)), "bam2freec_len requires a BAM file, got %s" % in_file
    base, _ = os.path.splitext(in_file)
    out_file = base + ".len"
    if file_exists(out_file):
        return out_file

    header = pysam.Samfile(in_file, 'rb').header
    with file_transaction(out_file) as tmp_out_file:
        with open(tmp_out_file, 'w') as out_handle:
            count = 1
            for line in header['SQ']:
                out_handle.write("\t".join([str(count), line['SN'],
                                            str(line['LN'])]) + "\n")
                count += 1
    return out_file

def is_sam(in_file):
    _, ext = os.path.splitext(in_file)
    if ext == ".sam":
        return True
    else:
        return False

def is_bam(in_file):
    _, ext = os.path.splitext(in_file)
    if ext == ".bam":
        return True
    else:
        return False

def make_refflat(genome_dir):
    """
    makes a refflat file for use with Picard from a GTF file
    """
    gtf_file = get_transcript_gtf(genome_dir)
    base, _ = os.path.splitext(gtf_file)
    refflat_file = base + ".refFlat"
    print "Making %s into a refFlat file named %s." % (gtf_file, refflat_file)
    if file_exists(refflat_file):
        print "%s already exists, skipping." % refflat_file
        return refflat_file

    with tmpfile(dir=os.getcwd(), prefix="genepred") as tmp_file:
        cmd = "gtfToGenePred {gtf_file} {tmp_file}".format(**locals())
        subprocess.check_call(cmd, shell=True)
        with open(tmp_file) as tmp_handle, open(refflat_file, "w") as out_handle:
            for line in tmp_handle:
                l = line.split("\t")
                l = [l[0]] + l
                out_handle.write("\t".join(l) + "\n")
    return refflat_file
