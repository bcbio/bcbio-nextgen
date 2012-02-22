#!/usr/bin/env python
"""Convert Hydra BEDPE output into VCF 4.1 format.

File format definitions:

http://code.google.com/p/hydra-sv/wiki/FileFormats
http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/
vcf-variant-call-format-version-41

Requirements:

PyVCF: https://github.com/jamescasbon/PyVCF
bx-python: https://bitbucket.org/james_taylor/bx-python/wiki/Home

Usage:
  hydra_to_vcf.py <hydra output file> <Genome in UCSC 2bit format>
"""
import os
import csv
import sys
import unittest
from collections import namedtuple

from bx.seq import twobit
import vcf

def main(hydra_file, genome_file):
    out_file = "{0}.vcf".format(os.path.splitext(hydra_file)[0])
    genome_2bit = twobit.TwoBitFile(open(genome_file))

def _vcf_info(start, end, mate_id):
    """Return breakend information line with mate and imprecise location.
    """
    return "SVTYPE=BND;MATEID={mate};IMPRECISE;CIPOS=0,{size}".format(
        mate=mate_id, size=end-start)

def _vcf_alt(base, other_chr, other_pos, other_isrc, is_first):
    if is_first:
        pipe = "]" if other_isrc else "["
        out_str = "{base}{pipe}{chr}:{pos}{pipe}"
    else:
        pipe = "[" if other_isrc else "]"
        out_str = "{pipe}{chr}:{pos}{pipe}{base}"
    return out_str.format(pipe=pipe, chr=other_chr, pos=other_pos,
                          base=base)

def build_vcf_parts(feature, genome_2bit):
    """Convert BedPe feature information into VCF part representation.

    Each feature will have two VCF lines for each side of the breakpoint.
    """
    VcfBreakend = namedtuple('VcfBreakend', ["chrom", "pos", "id", "ref", "alt",
                                             "info"])
    base1 = genome_2bit[feature.chrom1].get(
        feature.start1, feature.start1 + 1).upper()
    id1 = "hydra{0}a".format(feature.name)
    base2 = genome_2bit[feature.chrom2].get(
        feature.start2, feature.start2 + 1).upper()
    id2 = "hydra{0}b".format(feature.name)
    return (VcfBreakend(feature.chrom1, feature.start1, id1, base1,
                        _vcf_alt(base1, feature.chrom2, feature.start2,
                                 feature.strand2 == "-", True),
                        _vcf_info(feature.start1, feature.end1, id2)),
            VcfBreakend(feature.chrom2, feature.start2, id2, base2,
                        _vcf_alt(base2, feature.chrom1, feature.start1,
                                 feature.strand1 == "-", False),
                        _vcf_info(feature.start2, feature.end2, id1)))

def hydra_parser(in_file):
    """Parse hydra input file into namedtuple of values.
    """
    BedPe = namedtuple('BedPe', ["chrom1", "start1", "end1",
                                 "chrom2", "start2", "end2",
                                 "name", "strand1", "strand2"])
    with open(in_file) as in_handle:
        reader = csv.reader(in_handle, dialect="excel-tab")
        for line in reader:
            yield BedPe(line[0], int(line[1]), int(line[2]),
                        line[3], int(line[4]), int(line[5]),
                        line[6], line[8], line[9])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "Incorrect arguments"
        print __doc__
        sys.exist()
    main(sys.argv[1:])

# ## Test code

class HydraConvertTest(unittest.TestCase):
    """Test Hydra output conversion to VCF.
    """
    def setUp(self):
        self.work_dir = os.path.join(os.path.dirname(__file__), os.pardir,
                                     "tests", "data")
        self.in_file = os.path.join(self.work_dir, "structural",
                                    "NA12878-hydra.txt")
        self.genome_file = os.path.join(self.work_dir, "genomes", "hg19",
                                        "ucsc", "hg19.2bit")

    def test_1_input_parser(self):
        """Parse input file as BEDPE.
        """
        breakend = hydra_parser(self.in_file).next()
        assert breakend.chrom1 == "chr22"
        assert breakend.start1 == 9763 
        assert breakend.strand2 == "-"
        assert breakend.name == "2"

    def test_2_vcf_parts(self):
        """Convert BEDPE input line into VCF output parts.
        """
        genome_2bit = twobit.TwoBitFile(open(self.genome_file))
        breakend = hydra_parser(self.in_file).next()
        brend1, brend2 = build_vcf_parts(breakend, genome_2bit)
        assert brend1.alt == "G]chr22:10111]"
        assert brend2.info == "SVTYPE=BND;MATEID=hydra2a;IMPRECISE;CIPOS=0,102"
