"""Provide trimming of input reads from Fastq or BAM files.
"""
import os
import contextlib

from bcbio.utils import file_exists, save_diskspace, safe_makedir

from Bio.SeqIO.QualityIO import FastqGeneralIterator

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
    if (os.path.commonprefix([in_file, out_file]) ==
        os.path.split(os.path.dirname(out_file))[0]):
        utils.save_diskspace(in_file, "Trimmed to {}".format(out_file), config)

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
