"""
file format conversion utility functions
"""
import pysam
from bcbio.utils import (file_exists, replace_suffix)
from bcbio.distributed.transaction import file_transaction
import os

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

@expects("bam")
def bam2sam(in_file):
    """
    converts a bam file to a sam file
    bam2sam("file.bam") -> "file.sam"
    """
    out_file = replace_suffix(in_file, ".sam")
    if file_exists(out_file):
        return out_file
    with file_transaction(out_file) as tmp_out_file:
        pysam.view("-h", "-o" + tmp_out_file, in_file)
    return out_file

@expects("bam")
def bam2sizes(in_file):
    """
    converts a bam file to a chromosome sizes file
    chromosome sizes has the format:
    CHROM	SIZE

    For example:
    chr1	249250621
    chr2	243199373
    """
    base, _ = os.path.splitext(in_file)
    out_file = base + ".sizes"
    if file_exists(out_file):
        return out_file

    header = pysam.Samfile(in_file, 'rb').header
    with open(out_file, 'w') as out_handle:
        for line in header['SQ']:
            out_handle.write("\t".join([line['SN'], str(line['LN'])]) + "\n")
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
