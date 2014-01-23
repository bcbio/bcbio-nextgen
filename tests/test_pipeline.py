"""Test individual components of the analysis pipeline.
"""
import os
import shutil
import subprocess
import unittest

from nose.plugins.attrib import attr

from bcbio import utils
from bcbio.bam import fastq
from bcbio.distributed.messaging import parallel_runner
from bcbio.pipeline.config_utils import load_config
from bcbio.provenance import programs
from bcbio.variation import vcfutils

class RunInfoTest(unittest.TestCase):
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), "data")

    @attr(speed=1)
    def test_programs(self):
        """Identify programs and versions used in analysis.
        """
        config = load_config(os.path.join(self.data_dir, "automated",
                                          "post_process-sample.yaml"))
        print programs._get_versions(config)

class VCFUtilTest(unittest.TestCase):
    """Test various utilities for dealing with VCF files.
    """
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), "data")
        self.var_dir = os.path.join(self.data_dir, "variants")
        self.combo_file = os.path.join(self.var_dir, "S1_S2-combined.vcf.gz")

    @attr(speed=1)
    @attr(combo=True)
    def test_1_parallel_vcf_combine(self):
        """Parallel combination of VCF files, split by chromosome.
        """
        files = [os.path.join(self.var_dir, "S1-variants.vcf"), os.path.join(self.var_dir, "S2-variants.vcf")]
        ref_file = os.path.join(self.data_dir, "genomes", "hg19", "seq", "hg19.fa")
        config = load_config(os.path.join(self.data_dir, "automated",
                                          "post_process-sample.yaml"))
        run_parallel = parallel_runner({"type": "local", "cores": 1}, {}, config)
        region_dir = os.path.join(self.var_dir, "S1_S2-combined-regions")
        if os.path.exists(region_dir):
            shutil.rmtree(region_dir)
        if os.path.exists(self.combo_file):
            os.remove(self.combo_file)
        vcfutils.parallel_combine_variants(files, self.combo_file, ref_file, config, run_parallel)
        for fname in files:
            if os.path.exists(fname + ".gz"):
                subprocess.check_call(["gunzip", fname + ".gz"])
            if os.path.exists(fname + ".gz.tbi"):
                os.remove(fname + ".gz.tbi")

    @attr(speed=1)
    @attr(combo=True)
    def test_2_vcf_exclusion(self):
        """Exclude samples from VCF files.
        """
        ref_file = os.path.join(self.data_dir, "genomes", "hg19", "seq", "hg19.fa")
        config = load_config(os.path.join(self.data_dir, "automated",
                                          "post_process-sample.yaml"))
        out_file = utils.append_stem(self.combo_file, "-exclude")
        to_exclude = ["S1"]
        if os.path.exists(out_file):
            os.remove(out_file)
        vcfutils.exclude_samples(self.combo_file, out_file, to_exclude, ref_file, config)

    @attr(speed=1)
    @attr(template=True)
    def test_3_find_fastq_pairs(self):
        """Ensure we can correctly find paired fastq files.
        """
        test_pairs = ["/path/to/input/D1HJVACXX_2_AAGAGATC_1.fastq",
                      "/path/to/input/D1HJVACXX_3_AAGAGATC_1.fastq",
                      "/path/2/input/D1HJVACXX_2_AAGAGATC_2.fastq",
                      "/path/2/input/D1HJVACXX_3_AAGAGATC_2.fastq"]
        out = fastq.combine_pairs(test_pairs)
        assert out[0] == ["/path/to/input/D1HJVACXX_2_AAGAGATC_1.fastq",
                          "/path/2/input/D1HJVACXX_2_AAGAGATC_2.fastq"], out[0]
        assert out[1] == ["/path/to/input/D1HJVACXX_3_AAGAGATC_1.fastq",
                          "/path/2/input/D1HJVACXX_3_AAGAGATC_2.fastq"], out[1]


        test_pairs = ["/path/to/input/Tester_1_fastq.txt",
                      "/path/to/input/Tester_2_fastq.txt"]
        out = fastq.combine_pairs(test_pairs)
        assert out[0] == test_pairs, out[0]
