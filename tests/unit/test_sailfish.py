from bcbio.rnaseq import sailfish
from bcbio.utils import file_exists
import os

TEST_RESOURCES = {"fq1": "../data/110907_ERP000591/1_110907_ERP000591_1_fastq.txt",
                  "fq2": "../data/110907_ERP000591/1_110907_ERP000591_2_fastq.txt",
                  "gtf": "../data/genomes/mm9/rnaseq/ref-transcripts.gtf",
                  "fasta": "../data/genomes/mm9/seq/mm9.fa"}

from nose.plugins.attrib import attr

@attr("unit-broken")
def test_sailfish_index():
    out_dir = sailfish.sailfish_index(TEST_RESOURCES["gtf"], TEST_RESOURCES["fasta"])
    out_file = os.path.join(out_dir, "kmerEquivClasses.bin")
    assert file_exists(out_file)

@attr("unit-broken")
def test_sailfish():
    out_dir = sailfish.sailfish(TEST_RESOURCES["fq1"], TEST_RESOURCES["fq2"],
                                "yolo", TEST_RESOURCES["gtf"],
                                TEST_RESOURCES["fasta"], "unstranded")

