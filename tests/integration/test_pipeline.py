"""Test individual components of the analysis pipeline.
"""
import os
import sys
import shutil
import subprocess

from bcbio import utils
from bcbio.bam import fastq
from bcbio.distributed import prun
from bcbio.pipeline.config_utils import load_config

from tests.conftest import make_workdir, get_post_process_yaml

import pytest

class TestRunInfo(object):

    @pytest.mark.speed1
    def test_programs(self, data_dir):
        """Identify programs and versions used in analysis.
        """
        from bcbio.provenance import programs
        with make_workdir() as workdir:
            config = load_config(get_post_process_yaml(data_dir, workdir))
            print(programs._get_versions(config))


class TestVCFUtil(object):
    """Test various utilities for dealing with VCF files.
    """
    @property
    def data_dir(self):
        return os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            "data"
        )

    @property
    def automated_dir(self):
        return os.path.join(self.data_dir, "automated")

    @property
    def var_dir(self):
        return os.path.join(self.data_dir, "variants")

    @property
    def combo_file(self):
        return os.path.join(self.var_dir, "S1_S2-combined.vcf.gz")

    @property
    def ref_file(self):
        return os.path.join(self.data_dir, "genomes", "hg19", "seq", "hg19.fa")

    @pytest.mark.speed1
    @pytest.mark.combo
    def test_1_parallel_vcf_combine(self):
        """Parallel combination of VCF files, split by chromosome.
        """
        from bcbio.variation import vcfutils
        files = [
            os.path.join(self.var_dir, "S1-variants.vcf"),
            os.path.join(self.var_dir, "S2-variants.vcf")
        ]
        with make_workdir() as workdir:
            config = load_config(
                get_post_process_yaml(self.automated_dir, workdir))
            config["algorithm"] = {}
        region_dir = os.path.join(self.var_dir, "S1_S2-combined-regions")
        if os.path.exists(region_dir):
            shutil.rmtree(region_dir)
        if os.path.exists(self.combo_file):
            os.remove(self.combo_file)
        reqs = {"type": "local", "cores": 1}
        with prun.start(reqs, [[config]], config) as run_parallel:
            vcfutils.parallel_combine_variants(
                files, self.combo_file, self.ref_file, config, run_parallel)
        for fname in files:
            if os.path.exists(fname + ".gz"):
                subprocess.check_call(["gunzip", fname + ".gz"])
            if os.path.exists(fname + ".gz.tbi"):
                os.remove(fname + ".gz.tbi")

    @pytest.mark.speed1
    @pytest.mark.combo
    def test_2_vcf_exclusion(self):
        """Exclude samples from VCF files.
        """
        from bcbio.variation import vcfutils
        with make_workdir() as workdir:
            config = load_config(
                get_post_process_yaml(self.automated_dir, workdir))
            config["algorithm"] = {}
        out_file = utils.append_stem(self.combo_file, "-exclude")
        to_exclude = ["S1"]
        if os.path.exists(out_file):
            os.remove(out_file)
        vcfutils.exclude_samples(
            self.combo_file, out_file, to_exclude, self.ref_file, config)

    @pytest.mark.speed1
    @pytest.mark.combo
    def test_3_vcf_split_combine(self):
        """Split a VCF file into SNPs and indels, then combine back together.
        """
        from bcbio.variation import vcfutils
        with make_workdir() as workdir:
            config = load_config(get_post_process_yaml(
                self.automated_dir, workdir))
            config["algorithm"] = {}
        fname = os.path.join(self.var_dir, "S1-variants.vcf")
        snp_file, indel_file = vcfutils.split_snps_indels(
            fname, self.ref_file, config)
        merge_file = "%s-merge%s.gz" % os.path.splitext(fname)
        vcfutils.combine_variant_files(
            [snp_file, indel_file], merge_file, self.ref_file, config)
        for f in [snp_file, indel_file, merge_file]:
            self._remove_vcf(f)

    def _remove_vcf(self, f):
        for ext in ["", ".gz", ".gz.tbi", ".tbi"]:
            if os.path.exists(f + ext):
                os.remove(f + ext)

    @pytest.mark.speed1
    @pytest.mark.combo
    def test_4_vcf_sample_select(self, install_test_files, data_dir):
        """Select a sample from a VCF file.
        """
        from bcbio.variation import vcfutils
        fname = os.path.join(self.var_dir, "S1_S2-combined.vcf.gz")
        out_file = "%s-sampleselect%s" % utils.splitext_plus(fname)
        out_file = vcfutils.select_sample(fname, "S2", out_file, {})
        self._remove_vcf(out_file)

    @pytest.mark.speed1
    @pytest.mark.combo
    def test_5_find_fastq_pairs(self, install_test_files, data_dir):
        """Ensure we can correctly find paired fastq files.
        """
        test_pairs = ["/path/to/input/D1HJVACXX_2_AAGAGATC_1.fastq",
                      "/path/to/input/D1HJVACXX_5_AAGAGATC_1.fastq",
                      "/path/2/input/D1HJVACXX_2_AAGAGATC_2.fastq",
                      "/path/2/input/D1HJVACXX_5_AAGAGATC_2.fastq"]
        out = fastq.combine_pairs(test_pairs)
        assert out[0] == ["/path/to/input/D1HJVACXX_2_AAGAGATC_1.fastq",
                          "/path/2/input/D1HJVACXX_2_AAGAGATC_2.fastq"], out[0]
        assert out[1] == ["/path/to/input/D1HJVACXX_5_AAGAGATC_1.fastq",
                          "/path/2/input/D1HJVACXX_5_AAGAGATC_2.fastq"], out[1]

        test_pairs = ["/path/to/input/Tester_1_fastq.txt",
                      "/path/to/input/Tester_2_fastq.txt"]
        out = fastq.combine_pairs(test_pairs)
        assert out[0] == test_pairs, out[0]
