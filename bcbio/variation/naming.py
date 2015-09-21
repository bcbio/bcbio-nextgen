"""Fix chromosome naming incompatibilities for common issues, like hg19/GRCh37.

Fixes issues relating to chr1 versus 1 naming.

Uses Devon Ryan's great collection of contig mappings:

https://github.com/dpryan79/ChromosomeMappings
"""
import os
import requests

from bcbio import utils
from bcbio.bam import ref
from bcbio.distributed.transaction import file_transaction
from bcbio.variation import vcfutils

# ## Cached results
GMAP = {}

# read_mapping("https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh37_ensembl2UCSC.txt")
GMAP["hg19"] = {'GL000219.1': 'chrUn_gl000219', 'GL000192.1':
                'chr1_gl000192_random', 'GL000236.1': 'chrUn_gl000236', 'GL000211.1':
                'chrUn_gl000211', 'GL000234.1': 'chrUn_gl000234', '20': 'chr20', '21': 'chr21',
                '22': 'chr22', 'GL000196.1': 'chr8_gl000196_random', 'GL000213.1':
                'chrUn_gl000213', 'GL000205.1': 'chr17_gl000205_random', '4': 'chr4',
                'GL000222.1': 'chrUn_gl000222', 'GL000215.1': 'chrUn_gl000215', '8': 'chr8',
                'GL000232.1': 'chrUn_gl000232', 'GL000242.1': 'chrUn_gl000242', 'GL000244.1':
                'chrUn_gl000244', 'GL000223.1': 'chrUn_gl000223', 'GL000229.1':
                'chrUn_gl000229', 'GL000240.1': 'chrUn_gl000240', 'X': 'chrX', 'GL000202.1':
                'chr11_gl000202_random', 'GL000217.1': 'chrUn_gl000217', 'GL000200.1':
                'chr9_gl000200_random', 'GL000230.1': 'chrUn_gl000230', 'GL000206.1':
                'chr17_gl000206_random', 'HSCHR6_MHC_QBL': 'chr6_qbl_hap6', 'HSCHR6_MHC_MANN':
                'chr6_mann_hap4', 'GL000237.1': 'chrUn_gl000237', 'GL000204.1':
                'chr17_gl000204_random', 'GL000235.1': 'chrUn_gl000235', 'HSCHR6_MHC_APD':
                'chr6_apd_hap1', 'HSCHR6_MHC_COX': 'chr6_cox_hap2', '3': 'chr3', '7': 'chr7',
                'GL000233.1': 'chrUn_gl000233', 'GL000221.1': 'chrUn_gl000221', 'GL000220.1':
                'chrUn_gl000220', 'GL000245.1': 'chrUn_gl000245', 'GL000228.1':
                'chrUn_gl000228', 'GL000231.1': 'chrUn_gl000231', 'MT': 'chrM',
                'HSCHR6_MHC_SSTO': 'chr6_ssto_hap7', 'GL000238.1': 'chrUn_gl000238',
                'GL000195.1': 'chr7_gl000195_random', 'GL000249.1': 'chrUn_gl000249', '2':
                'chr2', '6': 'chr6', 'GL000247.1': 'chrUn_gl000247', 'GL000199.1':
                'chr9_gl000199_random', 'HSCHR6_MHC_DBB': 'chr6_dbb_hap3', 'GL000246.1':
                'chrUn_gl000246', 'GL000225.1': 'chrUn_gl000225', 'HSCHR4_1': 'chr4_ctg9_hap1',
                'GL000227.1': 'chrUn_gl000227', '11': 'chr11', '10': 'chr10', '13': 'chr13',
                '12': 'chr12', '15': 'chr15', '14': 'chr14', '17': 'chr17', '16': 'chr16', '19':
                'chr19', '18': 'chr18', 'GL000193.1': 'chr4_gl000193_random', 'GL000210.1':
                'chr21_gl000210_random', 'GL000239.1': 'chrUn_gl000239', 'GL000191.1':
                'chr1_gl000191_random', 'HSCHR17_1': 'chr17_ctg5_hap1', 'GL000194.1':
                'chr4_gl000194_random', 'GL000212.1': 'chrUn_gl000212', 'GL000248.1':
                'chrUn_gl000248', 'GL000197.1': 'chr8_gl000197_random', '1': 'chr1', '5':
                'chr5', 'GL000208.1': 'chr19_gl000208_random', '9': 'chr9', 'GL000214.1':
                'chrUn_gl000214', 'GL000224.1': 'chrUn_gl000224', 'GL000243.1':
                'chrUn_gl000243', 'HSCHR6_MHC_MCF': 'chr6_mcf_hap5', 'GL000209.1':
                'chr19_gl000209_random', 'GL000203.1': 'chr17_gl000203_random', 'GL000226.1':
                'chrUn_gl000226', 'GL000241.1': 'chrUn_gl000241', 'Y': 'chrY', 'GL000201.1':
                'chr9_gl000201_random', 'GL000198.1': 'chr9_gl000198_random', 'GL000216.1':
                'chrUn_gl000216', 'GL000218.1': 'chrUn_gl000218', 'GL000207.1':
                'chr18_gl000207_random'}

#read_mapping("https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh37_UCSC2ensembl.txt")
GMAP["GRCh37"] = {'chr19_gl000208_random': 'GL000208.1',
                  'chr21_gl000210_random': 'GL000210.1', 'chr6_apd_hap1': 'HSCHR6_MHC_APD',
                  'chr13': '13', 'chr12': '12', 'chr11': '11', 'chr10': '10', 'chr17': '17',
                  'chr16': '16', 'chr15': '15', 'chr14': '14', 'chr19': '19', 'chr18': '18',
                  'chr9_gl000198_random': 'GL000198.1', 'chrUn_gl000239': 'GL000239.1',
                  'chrUn_gl000238': 'GL000238.1', 'chrUn_gl000233': 'GL000233.1',
                  'chrUn_gl000232': 'GL000232.1', 'chrUn_gl000231': 'GL000231.1',
                  'chrUn_gl000230': 'GL000230.1', 'chrUn_gl000237': 'GL000237.1',
                  'chrUn_gl000236': 'GL000236.1', 'chrUn_gl000235': 'GL000235.1',
                  'chrUn_gl000234': 'GL000234.1', 'chr6_qbl_hap6': 'HSCHR6_MHC_QBL',
                  'chr11_gl000202_random': 'GL000202.1', 'chr17_gl000206_random': 'GL000206.1',
                  'chr6_cox_hap2': 'HSCHR6_MHC_COX', 'chr4_gl000193_random': 'GL000193.1',
                  'chrUn_gl000248': 'GL000248.1', 'chrUn_gl000249': 'GL000249.1',
                  'chrUn_gl000246': 'GL000246.1', 'chrUn_gl000247': 'GL000247.1',
                  'chrUn_gl000244': 'GL000244.1', 'chrUn_gl000245': 'GL000245.1',
                  'chrUn_gl000242': 'GL000242.1', 'chrUn_gl000243': 'GL000243.1',
                  'chrUn_gl000240': 'GL000240.1', 'chrUn_gl000241': 'GL000241.1',
                  'chr17_gl000204_random': 'GL000204.1', 'chr17_ctg5_hap1': 'HSCHR17_1',
                  'chr17_gl000205_random': 'GL000205.1', 'chr9_gl000199_random': 'GL000199.1',
                  'chr9_gl000201_random': 'GL000201.1', 'chr8': '8', 'chr6_ssto_hap7':
                  'HSCHR6_MHC_SSTO', 'chr8_gl000197_random': 'GL000197.1', 'chr6_dbb_hap3':
                  'HSCHR6_MHC_DBB', 'chr7_gl000195_random': 'GL000195.1', 'chr1_gl000191_random':
                  'GL000191.1', 'chr4_ctg9_hap1': 'HSCHR4_1', 'chr3': '3', 'chr2': '2', 'chr1':
                  '1', 'chr17_gl000203_random': 'GL000203.1', 'chrUn_gl000225': 'GL000225.1',
                  'chrY': 'Y', 'chrX': 'X', 'chr9_gl000200_random': 'GL000200.1', 'chr9': '9',
                  'chrM': 'MT', 'chr8_gl000196_random': 'GL000196.1', 'chr6_mann_hap4':
                  'HSCHR6_MHC_MANN', 'chrUn_gl000211': 'GL000211.1', 'chrUn_gl000213':
                  'GL000213.1', 'chrUn_gl000212': 'GL000212.1', 'chrUn_gl000215': 'GL000215.1',
                  'chrUn_gl000214': 'GL000214.1', 'chrUn_gl000217': 'GL000217.1',
                  'chrUn_gl000216': 'GL000216.1', 'chrUn_gl000219': 'GL000219.1',
                  'chrUn_gl000218': 'GL000218.1', 'chr19_gl000209_random': 'GL000209.1', 'chr22':
                  '22', 'chr20': '20', 'chr21': '21', 'chr6_mcf_hap5': 'HSCHR6_MHC_MCF', 'chr7':
                  '7', 'chr6': '6', 'chr5': '5', 'chr4': '4', 'chrUn_gl000228': 'GL000228.1',
                  'chrUn_gl000229': 'GL000229.1', 'chr1_gl000192_random': 'GL000192.1',
                  'chrUn_gl000224': 'GL000224.1', 'chr4_gl000194_random': 'GL000194.1',
                  'chrUn_gl000226': 'GL000226.1', 'chrUn_gl000227': 'GL000227.1',
                  'chrUn_gl000220': 'GL000220.1', 'chrUn_gl000221': 'GL000221.1',
                  'chrUn_gl000222': 'GL000222.1', 'chrUn_gl000223': 'GL000223.1',
                  'chr18_gl000207_random': 'GL000207.1'}

def handle_synonyms(in_file, ref_file, genome_build, work_dir, data):
    """Potentially handle remapping synonymous chromosome names between builds.

    Handles tab delimited file formats like BED and VCF where the contig
    is in the first column.
    """
    if genome_build in GMAP:
        mappings = GMAP[genome_build]
        contigs = set([c.name for c in ref.file_contigs(ref_file)])
        out_file = os.path.join(work_dir, "%s-fixed_contigs%s" % utils.splitext_plus(os.path.basename(in_file)))
        if not utils.file_exists(out_file):
            if out_file.endswith(".gz"):
                out_file = out_file.replace(".gz", "")
                needs_bgzip = True
            else:
                needs_bgzip = False
            checked_file = "%s.checked" % utils.splitext_plus(out_file)[0]
            if not _matches_contigs(in_file, contigs, checked_file):
                with file_transaction(data, out_file) as tx_out_file:
                    _write_newname_file(in_file, tx_out_file, mappings)
                if needs_bgzip:
                    out_file = vcfutils.bgzip_and_index(out_file, data["config"])
                return out_file
    return in_file

def _write_newname_file(in_file, out_file, mappings):
    """Re-write an input file with contigs matching the correct reference.
    """
    with open(in_file) as in_handle:
        with open(out_file, "w") as out_handle:
            for line in in_handle:
                if line.startswith("#"):
                    out_handle.write(line)
                else:
                    parts = line.split("\t")
                    new_contig = mappings.get(parts[0])
                    if new_contig:
                        parts[0] = new_contig
                        out_handle.write("\t".join(parts))

def _matches_contigs(in_file, contigs, checked_file):
    """Check if the contigs in the input file match the defined contigs in the reference genome.
    """
    tocheck_contigs = 2
    if utils.file_exists(checked_file):
        with open(checked_file) as in_handle:
            return in_handle.read().strip() == "match"
    else:
        with utils.open_gzipsafe(in_file) as in_handle:
            to_check = set([])
            for line in in_handle:
                if not line.startswith("#"):
                    to_check.add(line.split()[0])
                if len(to_check) >= tocheck_contigs:
                    break
        with open(checked_file, "w") as out_handle:
            if any([c not in contigs for c in to_check]):
                out_handle.write("different")
                return False
            else:
                out_handle.write("match")
                return True

# ## Retrieval of mappings

def read_mapping(url):
    mappings = {}
    for line in requests.get(url).text.split("\n"):
        parts = line.strip().split()
        if len(parts) == 2:
            first, second = parts
            mappings[str(first)] = str(second)
    return mappings