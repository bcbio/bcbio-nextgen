"""Test individual components of the analysis pipeline.
"""
import os
import sys
import shutil
import subprocess
import unittest

from nose.plugins.attrib import attr

from bcbio import utils
from bcbio.bam import fastq
from bcbio.distributed import prun
from bcbio.ngsalign import alignprep
from bcbio.pipeline.config_utils import load_config
from bcbio.provenance import programs
from bcbio.variation import vcfutils

sys.path.append(os.path.dirname(__file__))
from test_automated_analysis import get_post_process_yaml, make_workdir

class RunInfoTest(unittest.TestCase):
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), "data")
        self.automated_dir = os.path.join(os.path.dirname(__file__), "data", "automated")

    @attr(speed=1)
    @attr(blah=True)
    def test_programs(self):
        """Identify programs and versions used in analysis.
        """
        with make_workdir() as workdir:
            config = load_config(get_post_process_yaml(self.automated_dir, workdir))
            print programs._get_versions(config)

class MiscTest(unittest.TestCase):
    """Additional unit test cases to run regularly to confirm code logic.
    """
    @attr(speed=1)
    def test_align_split_size(self):
        """Checks on logic for estimating align split size.
        """
        assert alignprep._pick_align_split_size(10, 5, 20, 50) == 20
        assert alignprep._pick_align_split_size(250, 5, 20, 50) == 20
        assert alignprep._pick_align_split_size(500, 5, 20, 50) == 40
        assert alignprep._pick_align_split_size(750, 5, 20, 50) == 60

class VCFUtilTest(unittest.TestCase):
    """Test various utilities for dealing with VCF files.
    """
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), "data")
        self.var_dir = os.path.join(self.data_dir, "variants")
        self.combo_file = os.path.join(self.var_dir, "S1_S2-combined.vcf.gz")
        self.automated_dir = os.path.join(os.path.dirname(__file__), "data", "automated")

    @attr(speed=1)
    @attr(combo=True)
    def test_1_parallel_vcf_combine(self):
        """Parallel combination of VCF files, split by chromosome.
        """
        files = [os.path.join(self.var_dir, "S1-variants.vcf"), os.path.join(self.var_dir, "S2-variants.vcf")]
        ref_file = os.path.join(self.data_dir, "genomes", "hg19", "seq", "hg19.fa")
        with make_workdir() as workdir:
            config = load_config(get_post_process_yaml(self.automated_dir, workdir))
            config["algorithm"] = {}
        region_dir = os.path.join(self.var_dir, "S1_S2-combined-regions")
        if os.path.exists(region_dir):
            shutil.rmtree(region_dir)
        if os.path.exists(self.combo_file):
            os.remove(self.combo_file)
        with prun.start({"type": "local", "cores": 1}, [[config]], config) as run_parallel:
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
        with make_workdir() as workdir:
            config = load_config(get_post_process_yaml(self.automated_dir, workdir))
            config["algorithm"] = {}
        out_file = utils.append_stem(self.combo_file, "-exclude")
        to_exclude = ["S1"]
        if os.path.exists(out_file):
            os.remove(out_file)
        vcfutils.exclude_samples(self.combo_file, out_file, to_exclude, ref_file, config)

    @attr(speed=1)
    @attr(combo=True)
    def test_3_vcf_split_combine(self):
        """Split a VCF file into SNPs and indels, then combine back together.
        """
        with make_workdir() as workdir:
            config = load_config(get_post_process_yaml(self.automated_dir, workdir))
            config["algorithm"] = {}
        ref_file = os.path.join(self.data_dir, "genomes", "hg19", "seq", "hg19.fa")
        fname = os.path.join(self.var_dir, "S1-variants.vcf")
        snp_file, indel_file = vcfutils.split_snps_indels(fname, ref_file, config)
        merge_file = "%s-merge%s.gz" % os.path.splitext(fname)
        vcfutils.combine_variant_files([snp_file, indel_file], merge_file, ref_file,
                                       config)
        for f in [snp_file, indel_file, merge_file]:
            self._remove_vcf(f)

    def _remove_vcf(self, f):
        for ext in ["", ".gz", ".gz.tbi", ".tbi"]:
            if os.path.exists(f + ext):
                os.remove(f + ext)

    @attr(speed=1)
    @attr(combo=True)
    def test_4_vcf_sample_select(self):
        """Select a sample from a VCF file.
        """
        fname = os.path.join(self.var_dir, "S1_S2-combined.vcf.gz")
        out_file = "%s-sampleselect%s" % utils.splitext_plus(fname)
        out_file = vcfutils.select_sample(fname, "S2", out_file, {})
        self._remove_vcf(out_file)

    @attr(speed=1)
    @attr(template=True)
    def test_5_find_fastq_pairs(self):
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
