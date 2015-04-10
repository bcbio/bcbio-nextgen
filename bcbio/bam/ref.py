"""Manipulation functionality to deal with reference files.
"""
import collections

from bcbio import utils
from bcbio.pipeline import config_utils
from bcbio.provenance import do

def fasta_idx(in_file, config=None):
    """Retrieve samtools style fasta index.
    """
    fasta_index = in_file + ".fai"
    if not utils.file_exists(fasta_index):
        samtools = config_utils.get_program("samtools", config) if config else "samtools"
        cmd = "{samtools} faidx {in_file}"
        do.run(cmd.format(**locals()), "samtools faidx")
    return fasta_index

def file_contigs(ref_file, config=None):
    """Iterator of reference contigs and lengths from a reference file.
    """
    ContigInfo = collections.namedtuple("ContigInfo", "name size")
    with open(fasta_idx(ref_file, config)) as in_handle:
        for line in (l for l in in_handle if l.strip()):
            name, size = line.split()[:2]
            yield ContigInfo(name, int(size))
