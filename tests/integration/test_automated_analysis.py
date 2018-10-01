"""This directory is setup with configurations to run the main functional test.

It exercises a full analysis pipeline on a smaller subset of data.
"""
import glob
import os
import subprocess

import pytest
from tests.conftest import make_workdir
from tests.conftest import get_post_process_yaml


# Use 'install_required' mark to skip tests that cannot be run on Travis CI,
# because they require installation of additional dependencies (e.g. GATK)

@pytest.mark.speed3
@pytest.mark.skip(reason='Multiplexing not supporting in latest versions')
def test_3_full_pipeline(install_test_files, data_dir):
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
def test_4_empty_fastq(install_test_files, data_dir):
    """Handle analysis of empty fastq inputs from failed runs.
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, os.pardir, "110221_empty_FC12345AAXX"),
              os.path.join(data_dir, "run_info-empty.yaml")]
        subprocess.check_call(cl)


@pytest.mark.rnaseq
@pytest.mark.stranded
@pytest.mark.install_required
def test_2_stranded(install_test_files, data_dir):
    """Run an RNA-seq analysis with TopHat and generate gene-level counts.
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, os.pardir, "test_stranded"),
              os.path.join(data_dir, "run_info-stranded.yaml")]
        subprocess.check_call(cl)


@pytest.mark.rnaseq
@pytest.mark.tophat
@pytest.mark.rnaseq_vc
@pytest.mark.install_required
def test_2_rnaseq(install_test_files, data_dir):
    """Run an RNA-seq analysis with TopHat and generate gene-level counts.
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, os.pardir, "110907_ERP000591"),
              os.path.join(data_dir, "run_info-rnaseq.yaml")]
        subprocess.check_call(cl)


@pytest.mark.fusion
def test_2_fusion(install_test_files, data_dir):
    """Run an RNA-seq analysis and test fusion genes
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, os.pardir, "test_fusion"),
              os.path.join(data_dir, "run_info-fusion.yaml")]
        subprocess.check_call(cl)


@pytest.mark.star
@pytest.mark.rnaseq
@pytest.mark.rnaseq_standard
def test_2_star(install_test_files, data_dir):
    """Run an RNA-seq analysis with STAR and generate gene-level counts.
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, os.pardir, "test_fusion"),
              os.path.join(data_dir, "run_info-star.yaml")]
        subprocess.check_call(cl)


@pytest.mark.fastrnaseq
@pytest.mark.rnaseq
def test_2_fastrnaseq(install_test_files, data_dir):
    """Run a fast RNA-seq analysis
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, os.pardir, "test_fusion"),
              os.path.join(data_dir, "run_info-fastrnaseq.yaml")]
        subprocess.check_call(cl)


@pytest.mark.rnaseq
@pytest.mark.scrnaseq
def test_2_scrnaseq(install_test_files, data_dir):
    """Run a single-cell RNA-seq analysis
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, os.pardir, "Harvard-inDrop"),
              os.path.join(data_dir, "run_info-scrnaseq.yaml")]
        subprocess.check_call(cl)


@pytest.mark.rnaseq
@pytest.mark.rnaseq_standard
@pytest.mark.hisat2
def test_2_hisat2(install_test_files, data_dir):
    """Run an RNA-seq analysis with hisat2 and generate gene-level counts.
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, os.pardir, "110907_ERP000591"),
              os.path.join(data_dir, "run_info-hisat2.yaml")]
        subprocess.check_call(cl)


@pytest.mark.rnaseq
@pytest.mark.singleend
@pytest.mark.explant
def test_explant(install_test_files, data_dir):
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


@pytest.mark.srnaseq
@pytest.mark.srnaseq_star
def test_srnaseq_star(install_test_files, data_dir):
    """Run an sRNA-seq analysis.
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, os.pardir, "test_srnaseq"),
              os.path.join(data_dir, "run_info-srnaseq_star.yaml")]
        subprocess.check_call(cl)


@pytest.mark.srnaseq
@pytest.mark.srnaseq_bowtie
def test_srnaseq_bowtie(install_test_files, data_dir):
    """Run an sRNA-seq analysis.
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, os.pardir, "test_srnaseq"),
              os.path.join(data_dir, "run_info-srnaseq_bowtie.yaml")]
        subprocess.check_call(cl)


@pytest.mark.chipseq
def test_chipseq(install_test_files, data_dir):
    """
    Run a chip-seq alignment with Bowtie2
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, os.pardir, "test_chipseq"),
              os.path.join(data_dir, "run_info-chipseq.yaml")]
        subprocess.check_call(cl)


@pytest.mark.speed1
@pytest.mark.ensemble
@pytest.mark.install_required
def test_1_variantcall(install_test_files, data_dir):
    """Test variant calling with GATK pipeline.
    Requires GATK.
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, os.pardir, "100326_FC6107FAAXX"),
              os.path.join(data_dir, "run_info-variantcall.yaml")]
        subprocess.check_call(cl)


@pytest.mark.devel
@pytest.mark.speed1
def test_5_bam(install_test_files, data_dir):
    """Allow BAM files as input to pipeline.
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, "run_info-bam.yaml")]
        subprocess.check_call(cl)


@pytest.mark.speed2
@pytest.mark.install_required
def test_6_bamclean(install_test_files, data_dir):
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, os.pardir, "100326_FC6107FAAXX"),
              os.path.join(data_dir, "run_info-bamclean.yaml")]
        subprocess.check_call(cl)


@pytest.mark.speed2
@pytest.mark.cancer
@pytest.mark.cancermulti
@pytest.mark.install_required
def test_7_cancer(install_test_files, data_dir):
    """Test paired tumor-normal calling using multiple
    calling approaches: MuTect, VarScan, FreeBayes.
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, "run_info-cancer.yaml")]
        subprocess.check_call(cl)


@pytest.mark.cancer
@pytest.mark.cancerpanel
@pytest.mark.install_required
def test_7_cancer_nonormal(install_test_files, data_dir):
    """Test cancer calling without normal samples or with normal VCF panels.
    Requires MuTect and GATK.
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, "run_info-cancer2.yaml")]
        subprocess.check_call(cl)

@pytest.mark.cancer
@pytest.mark.cancerprecall
@pytest.mark.install_required
def test_7b_cancer_precall(install_test_files, data_dir):
    """Test somatic prioritization and effects prediction with pre-called inputs.
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, "run_info-cancer3.yaml")]
        subprocess.check_call(cl)

@pytest.mark.speed1
@pytest.mark.template
def test_8_template(install_test_files, data_dir):
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


@pytest.mark.joint
@pytest.mark.install_required
def test_9_joint(install_test_files, data_dir):
    """Perform joint calling/backfilling/squaring off following variant calling.
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, "run_info-joint.yaml")]
        subprocess.check_call(cl)

@pytest.mark.umibarcode
@pytest.mark.install_required
def test_10_umi(install_test_files, data_dir):
    """Allow BAM files as input to pipeline.
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, "run_info-umi.yaml")]
        subprocess.check_call(cl)

@pytest.mark.hla
def test_11_hla(install_test_files, data_dir):
    """Test HLA typing with OptiType.
    """
    from bcbio.hla import optitype
    hla_dir = os.path.join(data_dir, os.pardir, "100326_FC6107FAAXX", "hla")
    with make_workdir() as workdir:
        data = {"dirs": {"work": workdir},
                "rgnames": {"sample": "test"},
                "config": {},
                "hla": {"fastq": glob.glob(os.path.join(hla_dir, "*"))}}
        out = optitype.run(data)
        with open(out["hla"]["call_file"]) as in_handle:
            header = in_handle.readline().strip().split(",")
            hla_a = dict(zip(header), in_handle.readline().strip().split(","))
            assert hla_a["alleles"] == "HLA-A*11:01;HLA-A*24:02", hla_a
