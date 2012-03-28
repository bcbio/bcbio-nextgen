
#!/usr/bin/env python
"""Convert Hydra BEDPE output into VCF 4.1 format.

File format definitions:

http://code.google.com/p/hydra-sv/wiki/FileFormats
http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/
vcf-variant-call-format-version-41

Figure 2 of this review:

Quinlan, AR and IM Hall. 2011. Characterizing complex structural variation
in germline and somatic genomes. Trends in Genetics.
http://download.cell.com/trends/genetics/pdf/PIIS0168952511001685.pdf

Is a great overview of different structural variations and their breakpoint
representation.

Requirements:

bx-python: https://bitbucket.org/james_taylor/bx-python/wiki/Home

Usage:
  hydra_to_vcf.py <hydra output file> <Genome in UCSC 2bit format>
                  [--minsupport Minimum weighted support to include breakends]
"""
import os
import csv
import sys
import unittest
from collections import namedtuple
from operator import attrgetter
from optparse import OptionParser

from bx.seq import twobit

def main(hydra_file, genome_file, min_support=0):
    options = {"min_support": min_support}
    out_file = "{0}.vcf".format(os.path.splitext(hydra_file)[0])
    genome_2bit = twobit.TwoBitFile(open(genome_file))
    with open(out_file, "w") as out_handle:
        hydra_to_vcf_writer(hydra_file, genome_2bit, options, out_handle)

# ## Build VCF representation from Hydra BedPe format

def _vcf_info(start, end, mate_id):
    """Return breakend information line with mate and imprecise location.
    """
    return "SVTYPE=BND;MATEID={mate};IMPRECISE;CIPOS=0,{size}".format(
        mate=mate_id, size=end-start)

def _vcf_alt(base, other_chr, other_pos, isrc, is_first):
    """Create ALT allele line in VCF 4.1 format associating with other paired end.
    """
    if is_first:
        pipe = "[" if isrc else "]"
        out_str = "{base}{pipe}{chr}:{pos}{pipe}"
    else:
        pipe = "]" if isrc else "["
        out_str = "{pipe}{chr}:{pos}{pipe}{base}"
    return out_str.format(pipe=pipe, chr=other_chr, pos=other_pos + 1,
                          base=base)

def _breakend_orientation(strand1, strand2):
    """Convert BEDPE strand representation of breakpoints into VCF.

    | strand1  |  strand2 |     VCF      |
    +----------+----------+--------------+
    |   +      |     -    | t[p[ ]p]t    |
    |   +      |     +    | t]p] t]p]    |
    |   -      |     -    | [p[t [p[t    |
    |   -      |     +    | ]p]t t[p[    |
    """
    EndOrientation = namedtuple("EndOrientation",
                                ["is_first1", "is_rc1", "is_first2", "is_rc2"])
    if strand1 == "+" and strand2 == "-":
        return EndOrientation(True, True, False, True)
    elif strand1 == "+" and strand2 == "+":
        return EndOrientation(True, False, True, False)
    elif strand1 == "-" and strand2 == "-":
        return EndOrientation(False, False, False, False)
    elif strand1 == "-" and strand2 == "+":
        return EndOrientation(False, True, True, True)
    else:
        raise ValueError("Unexpected strand pairing: {0} {1}".format(
            strand1, strand2))

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
    orientation = _breakend_orientation(feature.strand1, feature.strand2)
    return (VcfBreakend(feature.chrom1, feature.start1, id1, base1,
                        _vcf_alt(base1, feature.chrom2, feature.start2,
                                 orientation.is_rc1, orientation.is_first1),
                        _vcf_info(feature.start1, feature.end1, id2)),
            VcfBreakend(feature.chrom2, feature.start2, id2, base2,
                        _vcf_alt(base2, feature.chrom1, feature.start1,
                                 orientation.is_rc2, orientation.is_first2),
                        _vcf_info(feature.start2, feature.end2, id1)))

# ## Parse Hydra output into BedPe tuple representation

def hydra_parser(in_file, options):
    """Parse hydra input file into namedtuple of values.
    """
    BedPe = namedtuple('BedPe', ["chrom1", "start1", "end1",
                                 "chrom2", "start2", "end2",
                                 "name", "strand1", "strand2",
                                 "support"])
    with open(in_file) as in_handle:
        reader = csv.reader(in_handle, dialect="excel-tab")
        for line in reader:
            cur = BedPe(line[0], int(line[1]), int(line[2]),
                        line[3], int(line[4]), int(line[5]),
                        line[6], line[8], line[9],
                        float(line[18]))
            if cur.support >= options["min_support"]:
                yield cur

# ## Write VCF output

def _write_vcf_header(out_handle):
    """Write VCF header information for Hydra structural variant.
    """
    def w(line):
        out_handle.write("{0}\n".format(line))
    w('##fileformat=VCFv4.1')
    w('##INFO=<ID=CIPOS,Number=2,Type=Integer,'
      'Description="Confidence interval around POS for imprecise variants">')
    w('##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">')
    w('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
    w('##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">')
    w('##source=hydra')
    w("#" + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]))

def _write_vcf_breakend(brend, out_handle):
    """Write out a single VCF line with breakpoint information.
    """
    out_handle.write("{0}\n".format("\t".join(str(x) for x in
        [brend.chrom, brend.pos + 1, brend.id, brend.ref, brend.alt,
         ".", "PASS", brend.info])))

def _get_vcf_breakends(hydra_file, genome_2bit, options):
    """Parse BEDPE input, yielding VCF ready breakends.
    """
    for feature in hydra_parser(hydra_file, options):
        for brend in build_vcf_parts(feature, genome_2bit):
            yield brend

def hydra_to_vcf_writer(hydra_file, genome_2bit, options, out_handle):
    """Write hydra output as sorted VCF file.

    Requires loading the hydra file into memory to perform sorting
    on output VCF. Could generalize this to no sorting or by-chromosome
    approach if this proves too memory intensive.
    """
    _write_vcf_header(out_handle)
    brends = list(_get_vcf_breakends(hydra_file, genome_2bit, options))
    brends.sort(key=attrgetter("chrom", "pos"))
    for brend in brends:
        _write_vcf_breakend(brend, out_handle)

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-s", "--minsupport", dest="minsupport", default=0)
    (options, args) = parser.parse_args()
    if len(args) != 2:
        print "Incorrect arguments"
        print __doc__
        sys.exist()
    main(args[0], args[1], int(options.minsupport))

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
        assert breakend.strand2 == "+"
        assert breakend.name == "2"
        assert breakend.support == 4.0

    def test_2_vcf_parts(self):
        """Convert BEDPE input line into VCF output parts.
        """
        genome_2bit = twobit.TwoBitFile(open(self.genome_file))
        breakends = hydra_parser(self.in_file)
        brend1, brend2 = build_vcf_parts(breakends.next(), genome_2bit)
        assert brend1.alt == "G]chr22:10112]"
        assert brend2.alt == "C]chr22:9764]"
        assert brend2.info == "SVTYPE=BND;MATEID=hydra2a;IMPRECISE;CIPOS=0,102"
        brend1, brend2 = build_vcf_parts(breakends.next(), genome_2bit)
        assert brend1.alt == "A[chr22:12112["
        assert brend2.alt == "]chr22:7764]G"
        brend1, brend2 = build_vcf_parts(breakends.next(), genome_2bit)
        assert brend1.alt == "[chr22:11112[A"
        assert brend2.alt == "[chr22:8764[T"
        brend1, brend2 = build_vcf_parts(breakends.next(), genome_2bit)
        assert brend1.alt == "]chr22:13112]G", brend1.alt
        assert brend2.alt == "A[chr22:9764[", brend2.alt
