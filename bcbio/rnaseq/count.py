"""
count number of reads mapping to features of transcripts

"""
import os
import sys
import HTSeq
import itertools
import pandas as pd

from bcbio.utils import (which, file_exists, get_in, safe_makedir)
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.log import logger
from bcbio import bam


def _get_files(data):
#    in_file = _get_sam_file(data)
    in_file = bam.sort(data["work_bam"], data["config"], order="queryname")
    gtf_file = data["genome_resources"]["rnaseq"]["transcripts"]
    work_dir = data["dirs"].get("work", "work")
    out_dir = os.path.join(work_dir, "htseq-count")
    out_file = os.path.join(out_dir, data['rgnames']['sample']) + ".counts"
    stats_file = os.path.join(out_dir, data['rgnames']['sample']) + ".stats"
    return in_file, gtf_file, out_file, stats_file


def is_countfile(in_file):
    with open(in_file) as in_handle:
        firstline = in_handle.next().split("\t")
    if len(firstline) != 2:
        return False
    try:
        int(firstline[1])
    except ValueError:
        return False
    return True


def _get_sam_file(data):
    in_file = data["work_bam"]
    config = data["config"]
    sorted = bam.sort(in_file, config, "queryname")
    sam = bam.bam_to_sam(sorted, config)
    return sam


def invert_strand(iv):
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError("Illegal strand")
    return iv2


class UnknownChrom(Exception):
    pass

def _get_stranded_flag(config):
    strand_flag = {"unstranded": "no",
                   "firststrand": "reverse",
                   "secondstrand": "yes"}
    stranded = _get_strandedness(config)
    assert stranded in strand_flag, ("%s is not a valid strandedness value. "
                                     "Valid values are 'firststrand', 'secondstrand', "
                                     "and 'unstranded")
    return strand_flag[stranded]

def _get_strandedness(config):
    return get_in(config, ("algorithm", "strandedness"), "unstranded").lower()


def htseq_count(data):
    """ adapted from Simon Anders htseq-count.py script
    http://www-huber.embl.de/users/anders/HTSeq/doc/count.html
    """

    sam_filename, gff_filename, out_file, stats_file = _get_files(data)
    stranded = _get_stranded_flag(data["config"])
    overlap_mode = "union"
    feature_type = "exon"
    id_attribute = "gene_id"
    minaqual = 0


    if file_exists(out_file):
        return out_file

    logger.info("Counting reads mapping to exons in %s using %s as the "
                    "annotation and strandedness as %s." % (os.path.basename(sam_filename),
                    os.path.basename(gff_filename), _get_strandedness(data["config"])))

    features = HTSeq.GenomicArrayOfSets("auto", stranded != "no")
    counts = {}

    # Try to open samfile to fail early in case it is not there
    open(sam_filename).close()

    gff = HTSeq.GFF_Reader(gff_filename)
    i = 0
    try:
        for f in gff:
            if f.type == feature_type:
                try:
                    feature_id = f.attr[id_attribute]
                except KeyError:
                    sys.exit("Feature %s does not contain a '%s' attribute" %
                             (f.name, id_attribute))
                if stranded != "no" and f.iv.strand == ".":
                    sys.exit("Feature %s at %s does not have strand "
                             "information but you are running htseq-count "
                             "in stranded mode. Use '--stranded=no'." %
                             (f.name, f.iv))
                features[f.iv] += feature_id
                counts[f.attr[id_attribute]] = 0
            i += 1
            if i % 100000 == 0:
                sys.stderr.write("%d GFF lines processed.\n" % i)
    except:
        sys.stderr.write("Error occured in %s.\n"
                         % gff.get_line_number_string())
        raise

    sys.stderr.write("%d GFF lines processed.\n" % i)

    if len(counts) == 0:
        sys.stderr.write("Warning: No features of type '%s' found.\n"
                         % feature_type)

    try:
        if bam.is_sam(sam_filename):
            read_seq = HTSeq.SAM_Reader(sam_filename)
        elif bam.is_bam(sam_filename):
            read_seq = HTSeq.BAM_Reader(sam_filename)
        first_read = iter(read_seq).next()
        pe_mode = first_read.paired_end
    except:
        sys.stderr.write("Error occured when reading first line of sam "
                         "file.\n")
        raise

    try:
        if pe_mode:
            read_seq_pe_file = read_seq
            read_seq = HTSeq.pair_SAM_alignments(read_seq)
        empty = 0
        ambiguous = 0
        notaligned = 0
        lowqual = 0
        nonunique = 0
        i = 0
        for r in read_seq:
            i += 1
            if not pe_mode:
                if not r.aligned:
                    notaligned += 1
                    continue
                try:
                    if r.optional_field("NH") > 1:
                        nonunique += 1
                        continue
                except KeyError:
                    pass
                if r.aQual < minaqual:
                    lowqual += 1
                    continue
                if stranded != "reverse":
                    iv_seq = (co.ref_iv for co in r.cigar if co.type == "M"
                              and co.size > 0)
                else:
                    iv_seq = (invert_strand(co.ref_iv) for co in r.cigar if
                              co.type == "M" and co.size > 0)
            else:
                if r[0] is not None and r[0].aligned:
                    if stranded != "reverse":
                        iv_seq = (co.ref_iv for co in r[0].cigar if
                                  co.type == "M" and co.size > 0)
                    else:
                        iv_seq = (invert_strand(co.ref_iv) for co in r[0].cigar if
                                  co.type == "M" and co.size > 0)
                else:
                    iv_seq = tuple()
                if r[1] is not None and r[1].aligned:
                    if stranded != "reverse":
                        iv_seq = itertools.chain(iv_seq,
                                                 (invert_strand(co.ref_iv) for co
                                                  in r[1].cigar if co.type == "M"
                                                  and co.size > 0))
                    else:
                        iv_seq = itertools.chain(iv_seq,
                                                 (co.ref_iv for co in r[1].cigar
                                                  if co.type == "M" and co.size
                                                  > 0))
                else:
                    if (r[0] is None) or not (r[0].aligned):
                        notaligned += 1
                        continue
                try:
                    if (r[0] is not None and r[0].optional_field("NH") > 1) or \
                       (r[1] is not None and r[1].optional_field("NH") > 1):
                        nonunique += 1
                        continue
                except KeyError:
                    pass
                if (r[0] and r[0].aQual < minaqual) or (r[1] and
                                                        r[1].aQual < minaqual):
                    lowqual += 1
                    continue

            try:
                if overlap_mode == "union":
                    fs = set()
                    for iv in iv_seq:
                        if iv.chrom not in features.chrom_vectors:
                            raise UnknownChrom
                        for iv2, fs2 in features[iv].steps():
                            fs = fs.union(fs2)
                elif (overlap_mode == "intersection-strict" or
                      overlap_mode == "intersection-nonempty"):
                    fs = None
                    for iv in iv_seq:
                        if iv.chrom not in features.chrom_vectors:
                            raise UnknownChrom
                        for iv2, fs2 in features[iv].steps():
                            if (len(fs2) > 0 or overlap_mode == "intersection-strict"):
                                if fs is None:
                                    fs = fs2.copy()
                                else:
                                    fs = fs.intersection(fs2)
                else:
                    sys.exit("Illegal overlap mode.")
                if fs is None or len(fs) == 0:
                    empty += 1
                elif len(fs) > 1:
                    ambiguous += 1
                else:
                    counts[list(fs)[0]] += 1
            except UnknownChrom:
                if not pe_mode:
                    rr = r
                else:
                    rr = r[0] if r[0] is not None else r[1]
                empty += 1

            if i % 100000 == 0:
                sys.stderr.write("%d sam %s processed.\n" % ( i, "lines " if not pe_mode else "line pairs"))

    except:
        if not pe_mode:
            sys.stderr.write("Error occured in %s.\n" % read_seq.get_line_number_string())
        else:
            sys.stderr.write("Error occured in %s.\n" % read_seq_pe_file.get_line_number_string() )
        raise

    sys.stderr.write("%d sam %s processed.\n" % (i, "lines " if not pe_mode else "line pairs"))

    with file_transaction(out_file) as tmp_out_file:
        with open(tmp_out_file, "w") as out_handle:
            on_feature = 0
            for fn in sorted(counts.keys()):
                on_feature += counts[fn]
                out_handle.write("%s\t%d\n" % (fn, counts[fn]))

    with file_transaction(stats_file) as tmp_stats_file:
        with open(tmp_stats_file, "w") as out_handle:
            out_handle.write("on_feature\t%d\n" % on_feature)
            out_handle.write("no_feature\t%d\n" % empty)
            out_handle.write("ambiguous\t%d\n" % ambiguous)
            out_handle.write("too_low_aQual\t%d\n" % lowqual)
            out_handle.write("not_aligned\t%d\n" % notaligned)
            out_handle.write("alignment_not_unique\t%d\n" % nonunique)

    return out_file

def combine_count_files(files, out_file=None):
    for f in files:
        assert file_exists(f), "%s does not exist or is empty."
        assert is_countfile(f), "%s does not seem to be a count file."
    col_names = [os.path.basename(os.path.splitext(x)[0]) for x in files]
    if not out_file:
        out_dir = os.path.join(os.path.dirname(files[0]))
        out_file = os.path.join(out_dir, "combined.counts")

    if file_exists(out_file):
        return out_file

    df = pd.io.parsers.read_table(f, sep="\t", index_col=0, header=None,
                                  names=[col_names[0]])
    for i, f in enumerate(files):
        if i == 0:
            df = pd.io.parsers.read_table(f, sep="\t", index_col=0, header=None,
                                          names=[col_names[0]])
        else:
            df = df.join(pd.io.parsers.read_table(f, sep="\t", index_col=0, header=None,
                                                  names=[col_names[i]]))

    df.to_csv(out_file, sep="\t", index_label="id")
    return out_file
