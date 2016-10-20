"""This directory is setup with configurations to run the main functional test.

It exercises a full analysis pipeline on a smaller subset of data.
"""
import os
import subprocess

import pytest
from tests.conftest import make_workdir
from tests.conftest import get_post_process_yaml


# Use 'install_required' mark to skip tests that cannot be run on Travis CI,
# because they require installation of additional dependencies (e.g. GATK)


class TestAutomatedAnalysis(object):
    """Setup a full automated analysis and run the pipeline.
    """

    @pytest.mark.skip(reason='Multiplexing not supporting in latest versions')
    @pytest.mark.speed3
    def test_3_full_pipeline(self, install_test_files, data_dir):
        """Run full automated analysis pipeline with multiplexing.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(data_dir, workdir),
                  os.path.join(data_dir, os.pardir, "110106_FC70BUKAAXX"),
                  os.path.join(data_dir, "run_info.yaml")]
            subprocess.check_call(cl)

    @pytest.mark.skip(reason='Multiplexing not supporting in latest versions')
    @pytest.mark.speed3
    def test_4_empty_fastq(self, install_test_files, data_dir):
        """Handle analysis of empty fastq inputs from failed runs.
        """
        with make_workdir() as workdir:
            cl = [
                "bcbio_nextgen.py",
                get_post_process_yaml(data_dir, workdir),
                os.path.join(data_dir, os.pardir, "110221_empty_FC12345AAXX"),
                os.path.join(data_dir, "run_info-empty.yaml")
            ]
            subprocess.check_call(cl)

    @pytest.marks('rnaseq', 'stranded', 'install_required')
    def test_2_stranded(self, install_test_files, data_dir):
        """Run an RNA-seq analysis with TopHat and generate gene-level counts.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(data_dir, workdir),
                  os.path.join(data_dir, os.pardir, "test_stranded"),
                  os.path.join(data_dir, "run_info-stranded.yaml")]
            subprocess.check_call(cl)

    @pytest.marks('rnaseq', 'tophat', 'install_required')
    def test_2_rnaseq(self, install_test_files, data_dir):
        """Run an RNA-seq analysis with TopHat and generate gene-level counts.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(data_dir, workdir),
                  os.path.join(data_dir, os.pardir, "110907_ERP000591"),
                  os.path.join(data_dir, "run_info-rnaseq.yaml")]
            subprocess.check_call(cl)

    @pytest.marks('fusion')
    def test_2_fusion(self, install_test_files, data_dir):
        """Run an RNA-seq analysis and test fusion genes
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(data_dir, workdir),
                  os.path.join(data_dir, os.pardir, "test_fusion"),
                  os.path.join(data_dir, "run_info-fusion.yaml")]
            subprocess.check_call(cl)

    @pytest.marks('star', 'rnaseq', 'rnaseq_stranded')
    def test_2_star(self, install_test_files, data_dir):
        """Run an RNA-seq analysis with STAR and generate gene-level counts.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(data_dir, workdir),
                  os.path.join(data_dir, os.pardir, "110907_ERP000591"),
                  os.path.join(data_dir, "run_info-star.yaml")]
            subprocess.check_call(cl)

    @pytest.marks('rnaseq', 'fastrnaseq')
    def test_2_fastrnaseq(self, install_test_files, data_dir):
        """Run a fast RNA-seq analysis
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(data_dir, workdir),
                  os.path.join(data_dir, os.pardir, "110907_ERP000591"),
                  os.path.join(data_dir, "run_info-fastrnaseq.yaml")]
            subprocess.check_call(cl)

    @pytest.marks('rnaseq', 'scrnaseq', 'install_required')
    def test_2_scrnaseq(self, install_test_files, data_dir):
        """Run a single-cell RNA-seq analysis
        Requires Umis.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(data_dir, workdir),
                  os.path.join(data_dir, os.pardir, "Harvard-inDrop"),
                  os.path.join(data_dir, "run_info-scrnaseq.yaml")]
            subprocess.check_call(cl)

    @pytest.marks('rnaseq', 'rnaseq_stranded', 'hisat2')
    def test_2_hisat2(self, install_test_files, data_dir):
        """Run an RNA-seq analysis with hisat2 and generate gene-level counts.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(data_dir, workdir),
                  os.path.join(data_dir, os.pardir, "110907_ERP000591"),
                  os.path.join(data_dir, "run_info-hisat2.yaml")]
            subprocess.check_call(cl)

    @pytest.marks('rnaseq', 'singleend', 'explant')
    def test_explant(self, install_test_files, data_dir):
        """
        Run an explant RNA-seq analysis with TopHat
        and generate gene-level counts.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(data_dir, workdir),
                  os.path.join(data_dir, os.pardir, "1_explant"),
                  os.path.join(data_dir, "run_info-explant.yaml")]
            subprocess.check_call(cl)

    @pytest.marks('srnaseq', 'srnaseq_star')
    def test_srnaseq_star(self, install_test_files, data_dir):
        """Run an sRNA-seq analysis.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(data_dir, workdir),
                  os.path.join(data_dir, os.pardir, "test_srnaseq"),
                  os.path.join(data_dir, "run_info-srnaseq_star.yaml")]
            subprocess.check_call(cl)

    @pytest.marks('srnaseq', 'srnaseq_bowtie')
    def test_srnaseq_bowtie(self, install_test_files, data_dir):
        """Run an sRNA-seq analysis.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(data_dir, workdir),
                  os.path.join(data_dir, os.pardir, "test_srnaseq"),
                  os.path.join(data_dir, "run_info-srnaseq_bowtie.yaml")]
            subprocess.check_call(cl)

    @pytest.marks('chipseq')
    def test_chipseq(self, install_test_files, data_dir):
        """
        Run a chip-seq alignment with Bowtie2
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(data_dir, workdir),
                  os.path.join(data_dir, os.pardir, "test_chipseq"),
                  os.path.join(data_dir, "run_info-chipseq.yaml")]
            subprocess.check_call(cl)

    @pytest.marks('speed1', 'ensemble', 'install_required')
    def test_1_variantcall(self, install_test_files, data_dir):
        """Test variant calling with GATK pipeline.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(data_dir, workdir),
                  os.path.join(data_dir, os.pardir, "100326_FC6107FAAXX"),
                  os.path.join(data_dir, "run_info-variantcall.yaml")]
            subprocess.check_call(cl)

    @pytest.marks('devel', 'speed1')
    def test_5_bam(self, install_test_files, data_dir):
        """Allow BAM files as input to pipeline.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(data_dir, workdir),
                  os.path.join(data_dir, "run_info-bam.yaml")]
            subprocess.check_call(cl)

    @pytest.marks('speed2', 'install_required')
    def test_6_bamclean(self, install_test_files, data_dir):
        """Clean problem BAM input files that do not require alignment.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(data_dir, workdir),
                  os.path.join(data_dir, os.pardir, "100326_FC6107FAAXX"),
                  os.path.join(data_dir, "run_info-bamclean.yaml")]
            subprocess.check_call(cl)

    @pytest.marks('speed2', 'cancer', 'cancermulti')
    def test_7_cancer(self, install_test_files, data_dir):
        """Test paired tumor-normal calling using multiple
        calling approaches: MuTect, VarScan, FreeBayes.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(data_dir, workdir),
                  os.path.join(data_dir, "run_info-cancer.yaml")]
            subprocess.check_call(cl)

    @pytest.marks('speed2', 'cancer', 'cancerpanel', 'install_required')
    def test_7_cancer_nonormal(self, install_test_files, data_dir):
        """Test cancer calling without normal samples or with normal VCF panels.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(data_dir, workdir),
                  os.path.join(data_dir, "run_info-cancer2.yaml")]
            subprocess.check_call(cl)

    @pytest.marks('speed1', 'template')
    def test_8_template(self, install_test_files, data_dir):
        """Create a project template from input files and metadata configuration.
        """
        fc_dir = os.path.join(data_dir, os.pardir, "100326_FC6107FAAXX")
        with make_workdir():
            cl = ["bcbio_nextgen.py", "-w", "template", "--only-metadata",
                  "freebayes-variant",
                  os.path.join(fc_dir, "100326.csv"),
                  os.path.join(fc_dir, "7_100326_FC6107FAAXX_1_fastq.txt"),
                  os.path.join(fc_dir, "7_100326_FC6107FAAXX_2_fastq.txt"),
                  os.path.join(fc_dir, "8_100326_FC6107FAAXX.bam")]
            subprocess.check_call(cl)

    @pytest.marks('joint', 'install_required')
    def test_9_joint(self, install_test_files, data_dir):
        """Perform joint calling/backfilling/squaring off following variant calling.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(data_dir, workdir),
                  os.path.join(data_dir, "run_info-joint.yaml")]
            subprocess.check_call(cl)
