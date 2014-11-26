import unittest
from nose.plugins.attrib import attr

import os
from bcbio.rnaseq import rsem
from bcbio.utils import file_exists


class Rsem(unittest.TestCase):

    def setUp(self):
        self.fasta = os.path.abspath("../data/genomes/mm9/seq/mm9.fa")
        self.gtf = os.path.abspath("../data/genomes/mm9/rnaseq/ref-transcripts.gtf")
        self.rsem_genome_dir = os.path.abspath("bcbio-nextgen-test-data/data/rsem")
        self.bam = os.path.abspath("bcbio-nextgen-test-data/data/rsem/Test1.transcriptome.bam")


    @attr("unit")
    def test_prepare_rsem_reference(self):
        reference = rsem.prepare_rsem_reference(self.gtf, self.fasta, "chr22")
        self.assertTrue(file_exists(os.path.join(reference, "chr22.chrlist")))

    @attr("unit")
    def test_rsem_calculate_expression(self):
        out_file = rsem.rsem_calculate_expression(self.bam, self.rsem_genome_dir, "Test",
                                                  "chr22", "tmp")
