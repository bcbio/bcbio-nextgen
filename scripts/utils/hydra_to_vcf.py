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
from collections import namedtuple, defaultdict
from operator import attrgetter
from optparse import OptionParser

from bx.seq import twobit
from bx.intervals.cluster import ClusterTree

def main(hydra_file, genome_file, min_support=0):
    options = {"min_support": min_support, "max_single_size": 10000}
    out_file = "{0}.vcf".format(os.path.splitext(hydra_file)[0])
    genome_2bit = twobit.TwoBitFile(open(genome_file))
    with open(out_file, "w") as out_handle:
        hydra_to_vcf_writer(hydra_file, genome_2bit, options, out_handle)

# ## Build VCF breakend representation from Hydra BedPe format

def _vcf_info(start, end, mate_id, info=None):
    """Return breakend information line with mate and imprecise location.
    """
    out = "SVTYPE=BND;MATEID={mate};IMPRECISE;CIPOS=0,{size}".format(
        mate=mate_id, size=end-start)
    if info is not None:
        extra_info = ";".join("{0}={1}".format(k, v) for k, v in info.iteritems())
        out = "{0};{1}".format(out, extra_info)
    return out

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

VcfLine = namedtuple('VcfLine', ["chrom", "pos", "id", "ref", "alt", "info"])

def build_vcf_parts(feature, genome_2bit, info=None):
    """Convert BedPe feature information into VCF part representation.

    Each feature will have two VCF lines for each side of the breakpoint.
    """
    base1 = genome_2bit[feature.chrom1].get(
        feature.start1, feature.start1 + 1).upper()
    id1 = "hydra{0}a".format(feature.name)
    base2 = genome_2bit[feature.chrom2].get(
        feature.start2, feature.start2 + 1).upper()
    id2 = "hydra{0}b".format(feature.name)
    orientation = _breakend_orientation(feature.strand1, feature.strand2)
    return (VcfLine(feature.chrom1, feature.start1, id1, base1,
                    _vcf_alt(base1, feature.chrom2, feature.start2,
                             orientation.is_rc1, orientation.is_first1),
                    _vcf_info(feature.start1, feature.end1, id2, info)),
            VcfLine(feature.chrom2, feature.start2, id2, base2,
                    _vcf_alt(base2, feature.chrom1, feature.start1,
                             orientation.is_rc2, orientation.is_first2),
                    _vcf_info(feature.start2, feature.end2, id1, info)))

# ## Represent standard variants types
# Convert breakends into deletions, tandem duplications and inversions

def is_deletion(x, options):
    strand_orientation = ["+", "-"]
    return (x.chrom1 == x.chrom2 and [x.strand1, x.strand2] == strand_orientation
            and (x.start2 - x.start1) < options.get("max_single_size", 0))

def _vcf_single_end_info(x, svtype, is_removal=False):
    if is_removal:
        length = x.start1 - x.start2
    else:
        length = x.start2 - x.start1
    return "SVTYPE={type};IMPRECISE;CIPOS=0,{size1};CIEND=0,{size2};" \
           "END={end};SVLEN={length}".format(size1=x.end1 - x.start1,
                                             size2=x.end2 - x.start2,
                                             end=x.start2,
                                             type=svtype,
                                             length=length)

def build_vcf_deletion(x, genome_2bit):
    """Provide representation of deletion from BedPE breakpoints.
    """
    base1 = genome_2bit[x.chrom1].get(x.start1, x.start1 + 1).upper()
    id1 = "hydra{0}".format(x.name)
    return VcfLine(x.chrom1, x.start1, id1, base1, "<DEL>",
                   _vcf_single_end_info(x, "DEL", True))

def is_tandem_dup(x, options):
    strand_orientation = ["-", "+"]
    return (x.chrom1 == x.chrom2 and [x.strand1, x.strand2] == strand_orientation
            and (x.start2 - x.start1) < options.get("max_single_size", 0))

def build_tandem_deletion(x, genome_2bit):
    """Provide representation of tandem duplication.
    """
    base1 = genome_2bit[x.chrom1].get(x.start1, x.start1 + 1).upper()
    id1 = "hydra{0}".format(x.name)
    return VcfLine(x.chrom1, x.start1, id1, base1, "<DUP:TANDEM>",
                   _vcf_single_end_info(x, "DUP"))

def is_inversion(x1, x2):
    strand1 = ["+", "+"]
    strand2 = ["-", "-"]
    return (x1.chrom1 == x1.chrom2 and x1.chrom1 == x2.chrom1 and
            (([x1.strand1, x1.strand2] == strand1 and
              [x2.strand1, x2.strand2] == strand2) or
             ([x1.strand1, x1.strand2] == strand2 and
              [x2.strand1, x2.strand2] == strand1)))

def build_vcf_inversion(x1, x2, genome_2bit):
    """Provide representation of inversion from BedPE breakpoints.
    """
    id1 = "hydra{0}".format(x1.name)
    start_coords = sorted([x1.start1, x1.end1, x2.start1, x2.end1])
    end_coords = sorted([x1.start2, x1.end2, x2.start2, x2.start2])
    start_pos = (start_coords[1] + start_coords[2]) // 2
    end_pos = (end_coords[1] + end_coords[2]) // 2
    base1 = genome_2bit[x1.chrom1].get(start_pos, start_pos + 1).upper()
    info = "SVTYPE=INV;IMPRECISE;CIPOS={cip1},{cip2};CIEND={cie1},{cie2};" \
           "END={end};SVLEN={length}".format(cip1=start_pos - start_coords[0],
                                             cip2=start_coords[-1] - start_pos,
                                             cie1=end_pos - end_coords[0],
                                             cie2=end_coords[-1] - end_pos,
                                             end=end_pos,
                                             length=end_pos-start_pos)
    return VcfLine(x1.chrom1, start_pos, id1, base1, "<INV>", info)

def is_translocation(x1, x2):
    strand1 = ["+", "-"]
    strand2 = ["-", "+"]
    return (x1.chrom1 != x1.chrom2 and
            ([x1.strand1, x1.strand2] == strand1 and
             [x2.strand1, x2.strand2] == strand2) or
            ([x1.strand1, x1.strand2] == strand2 and
             [x2.strand1, x2.strand2] == strand1))

def get_translocation_info(x1, x2):
    return {"EVENT": "translocation_{0}_{1}".format(x1.name, x2.name)}

# ## Parse Hydra output into BedPe tuple representation

def hydra_parser(in_file, options=None):
    """Parse hydra input file into namedtuple of values.
    """
    if options is None: options = {}
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
            if cur.support >= options.get("min_support", 0):
                yield cur

def _cluster_by(end_iter, attr1, attr2, cluster_distance):
    """Cluster breakends by specified attributes.
    """
    ClusterInfo = namedtuple("ClusterInfo", ["chroms", "clusters", "lookup"])
    chr_clusters = {}
    chroms = []
    brends_by_id = {}
    for brend in end_iter:
        if not chr_clusters.has_key(brend.chrom1):
            chroms.append(brend.chrom1)
            chr_clusters[brend.chrom1] = ClusterTree(cluster_distance, 1)
        brends_by_id[int(brend.name)] = brend
        chr_clusters[brend.chrom1].insert(getattr(brend, attr1),
                                          getattr(brend, attr2),
                                          int(brend.name))
    return ClusterInfo(chroms, chr_clusters, brends_by_id)

def _calculate_cluster_distance(end_iter):
    """Compute allowed distance for clustering based on end confidence intervals.
    """
    out = []
    sizes = []
    for x in end_iter:
        out.append(x)
        sizes.append(x.end1 - x.start1)
        sizes.append(x.end2 - x.start2)
    distance = sum(sizes) // len(sizes)
    return distance, out

def group_hydra_breakends(end_iter):
    """Group together hydra breakends with overlapping ends.

    This provides a way to identify inversions, translocations
    and insertions present in hydra break point ends. We cluster together the
    endpoints and return together any items with closely oriented pairs.
    This helps in describing more complex rearrangement events.
    """
    cluster_distance, all_ends = _calculate_cluster_distance(end_iter)
    first_cluster = _cluster_by(all_ends, "start1", "end1", cluster_distance)
    for chrom in first_cluster.chroms:
        for _, _, brends in first_cluster.clusters[chrom].getregions():
            if len(brends) == 1:
                yield [first_cluster.lookup[brends[0]]]
            else:
                second_cluster = _cluster_by([first_cluster.lookup[x] for x in brends],
                                             "start2", "end2", cluster_distance)
                for chrom2 in second_cluster.chroms:
                    for _, _, brends in second_cluster.clusters[chrom].getregions():
                        yield [second_cluster.lookup[x] for x in brends]

# ## Write VCF output

def _write_vcf_header(out_handle):
    """Write VCF header information for Hydra structural variant.
    """
    def w(line):
        out_handle.write("{0}\n".format(line))
    w('##fileformat=VCFv4.1')
    w('##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">')
    w('##INFO=<ID=END,Number=1,Type=Integer,'
      'Description="End position of the variant described in this record">')
    w('##INFO=<ID=CIPOS,Number=2,Type=Integer,'
      'Description="Confidence interval around POS for imprecise variants">')
    w('##INFO=<ID=CIEND,Number=2,Type=Integer,'
      'Description="Confidence interval around END for imprecise variants">')
    w('##INFO=<ID=SVLEN,Number=.,Type=Integer,'
      'Description="Difference in length between REF and ALT alleles">')
    w('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
    w('##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">')
    w('##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">')
    w('##ALT=<ID=DEL,Description="Deletion">')
    w('##ALT=<ID=INV,Description="Inversion">')
    w('##ALT=<ID=DUP,Description="Duplication">')
    w('##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">')
    w('##source=hydra')
    w("#" + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]))

def _write_vcf_breakend(brend, out_handle):
    """Write out a single VCF line with breakpoint information.
    """
    out_handle.write("{0}\n".format("\t".join(str(x) for x in
        [brend.chrom, brend.pos + 1, brend.id, brend.ref, brend.alt,
         ".", "PASS", brend.info])))

def _get_vcf_breakends(hydra_file, genome_2bit, options=None):
    """Parse BEDPE input, yielding VCF ready breakends.
    """
    if options is None: options = {}
    for features in group_hydra_breakends(hydra_parser(hydra_file, options)):
        if len(features) == 1 and is_deletion(features[0], options):
            yield build_vcf_deletion(features[0], genome_2bit)
        elif len(features) == 1 and is_tandem_dup(features[0], options):
            yield build_tandem_deletion(features[0], genome_2bit)
        elif len(features) == 2 and is_inversion(*features):
            yield build_vcf_inversion(features[0], features[1], genome_2bit)
        elif len(features) == 2 and is_translocation(*features):
            info = get_translocation_info(features[0], features[1])
            for feature in features:
                for brend in build_vcf_parts(feature, genome_2bit, info):
                    yield brend
        else:
            for feature in features:
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
        breakend = next(hydra_parser(self.in_file))
        assert breakend.chrom1 == "chr22"
        assert breakend.start1 == 9763 
        assert breakend.strand2 == "+"
        assert breakend.name == "1"
        assert breakend.support == 4.0

    def test_2_vcf_parts(self):
        """Convert BEDPE input line into VCF output parts.
        """
        genome_2bit = twobit.TwoBitFile(open(self.genome_file))
        breakends = hydra_parser(self.in_file)
        brend1, brend2 = build_vcf_parts(next(breakends), genome_2bit)
        assert brend1.alt == "G]chr22:10112]"
        assert brend2.alt == "C]chr22:9764]"
        assert brend2.info == "SVTYPE=BND;MATEID=hydra1a;IMPRECISE;CIPOS=0,102", brend2.info
        brend1, brend2 = build_vcf_parts(next(breakends), genome_2bit)
        assert brend1.alt == "A[chr22:12112["
        assert brend2.alt == "]chr22:7764]G"
        brend1, brend2 = build_vcf_parts(next(breakends), genome_2bit)
        assert brend1.alt == "[chr22:11112[A"
        assert brend2.alt == "[chr22:8764[T"
        brend1, brend2 = build_vcf_parts(next(breakends), genome_2bit)
        assert brend1.alt == "]chr22:13112]G", brend1.alt
        assert brend2.alt == "A[chr22:9764[", brend2.alt

    def test_3_deletions(self):
        """Convert BEDPE breakends that form a deletion.
        """
        genome_2bit = twobit.TwoBitFile(open(self.genome_file))
        parts = _get_vcf_breakends(self.in_file, genome_2bit, {"max_single_size": 5000})
        deletion = next(parts)
        assert deletion.alt == "<DEL>", deletion
        assert "SVLEN=-4348" in deletion.info
