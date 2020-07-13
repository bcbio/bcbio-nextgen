"""This directory is setup with configurations to run the main functional test.

It exercises a full analysis pipeline on a smaller subset of data.

Use 'install_required' mark to skip tests that cannot be run on Travis CI,
because they require installation of additional dependencies (e.g. GATK)
"""
import glob
import os
import subprocess

import pytest


@pytest.mark.speed3
@pytest.mark.xfail(reason='Multiplexing not supporting in latest versions', run=False)
def test_full_pipeline_with_multiplexing(install_test_files, data_dir, global_config):
    fc_dir = os.path.join(data_dir, os.pardir, '110106_FC70BUKAAXX')
    run_config = os.path.join(data_dir, 'run_info.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])


@pytest.mark.speed3
@pytest.mark.xfail(reason='Multiplexing not supporting in latest versions', run=False)
def test_handle_empty_fastq_inputs_from_failed_runs(install_test_files, data_dir, global_config):
    fc_dir = os.path.join(data_dir, os.pardir, '110221_empty_FC12345AAXX')
    run_config = os.path.join(data_dir, 'run_info-empty.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])


@pytest.mark.speed1
@pytest.mark.rnaseq
@pytest.mark.stranded
@pytest.mark.install_required
def test_stranded(install_test_files, data_dir, global_config):
    """Run an RNA-seq analysis with TopHat and generate gene-level counts"""
    fc_dir = os.path.join(data_dir, os.pardir, 'test_stranded')
    run_config = os.path.join(data_dir, 'run_info-stranded.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])


@pytest.mark.rnaseq
@pytest.mark.tophat
@pytest.mark.rnaseq_vc
@pytest.mark.install_required
@pytest.mark.skip(reason="tophat is no longer supported.")
def test_rnaseq_with_tophat(install_test_files, data_dir, global_config):
    fc_dir = os.path.join(data_dir, os.pardir, '110907_ERP000591')
    run_config = os.path.join(data_dir, 'run_info-rnaseq.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])


@pytest.mark.fusion
def test_rnaseq_fusion_genes(install_test_files, data_dir, global_config):
    fc_dir = os.path.join(data_dir, os.pardir, 'test_fusion')
    run_config = os.path.join(data_dir, 'run_info-fusion.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])


@pytest.mark.star
@pytest.mark.rnaseq
@pytest.mark.rnaseq_standard
def test_rnaseq_with_star(install_test_files, data_dir, global_config):
    fc_dir = os.path.join(data_dir, os.pardir, 'test_fusion')
    run_config = os.path.join(data_dir, 'run_info-star.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])


@pytest.mark.fastrnaseq
@pytest.mark.rnaseq
def test_fast_rnaseq(install_test_files, data_dir, global_config):
    fc_dir = os.path.join(data_dir, os.pardir, 'test_fusion')
    run_config = os.path.join(data_dir, 'run_info-fastrnaseq.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])


@pytest.mark.rnaseq
@pytest.mark.scrnaseq
def test_scrnaseq(install_test_files, data_dir, global_config):
    fc_dir = os.path.join(data_dir, os.pardir, 'Harvard-inDrop')
    run_config = os.path.join(data_dir, 'run_info-scrnaseq.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])


@pytest.mark.rnaseq
@pytest.mark.rnaseq_standard
@pytest.mark.hisat2
def test_rnaseq_with_hisat2(install_test_files, data_dir, global_config):
    fc_dir = os.path.join(data_dir, os.pardir, '110907_ERP000591')
    run_config = os.path.join(data_dir, 'run_info-hisat2.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])


@pytest.mark.rnaseq
@pytest.mark.singleend
@pytest.mark.explant
def test_explant_with_tophat(install_test_files, data_dir, global_config):
    fc_dir = os.path.join(data_dir, os.pardir, '1_explant')
    run_config = os.path.join(data_dir, 'run_info-explant.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])


@pytest.mark.srnaseq
@pytest.mark.srnaseq_star
def test_srnaseq_star(install_test_files, data_dir, global_config):
    fc_dir = os.path.join(data_dir, os.pardir, 'test_srnaseq')
    run_config = os.path.join(data_dir, 'run_info-srnaseq_star.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])


@pytest.mark.srnaseq
@pytest.mark.srnaseq_bowtie
def test_srnaseq_bowtie(install_test_files, data_dir, global_config):
    fc_dir = os.path.join(data_dir, os.pardir, 'test_srnaseq')
    run_config = os.path.join(data_dir, 'run_info-srnaseq_bowtie.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])


@pytest.mark.chipseq
@pytest.mark.xfail(reason='https://github.com/bcbio/bcbio-nextgen/issues/3224', run=False)
def test_chipseq_with_bowtie2(install_test_files, data_dir, global_config):
    fc_dir = os.path.join(data_dir, os.pardir, 'test_chipseq')
    run_config = os.path.join(data_dir, 'run_info-chipseq.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])


@pytest.mark.atacseq
@pytest.mark.xfail(reason='https://github.com/bcbio/bcbio-nextgen/issues/3225', run=False)
def test_atacseq(install_test_files, data_dir, global_config):
    fc_dir = os.path.join(data_dir, os.pardir, 'test_atacseq')
    run_config = os.path.join(data_dir, 'run_info-atacseq.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])


@pytest.mark.speed1
@pytest.mark.ensemble
@pytest.mark.install_required
def test_variant_call_gatk(install_test_files, global_config, data_dir):
    # requires GATK
    fc_dir = os.path.join(data_dir, os.pardir, '100326_FC6107FAAXX')
    run_config = os.path.join(data_dir, 'run_info-variantcall.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])


@pytest.mark.devel
@pytest.mark.speed1
def test_variant2_pipeline_with_bam_input(install_test_files, global_config, data_dir):
    run_config = os.path.join(data_dir, 'run_info-bam.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, run_config])


@pytest.mark.speed2
@pytest.mark.install_required
def test_bamclean(install_test_files, data_dir, global_config):
    fc_dir = os.path.join(data_dir, os.pardir, '100326_FC6107FAAXX')
    run_config = os.path.join(data_dir, 'run_info-bamclean.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])


@pytest.mark.speed2
@pytest.mark.cancer
@pytest.mark.cancermulti
@pytest.mark.install_required
def test_paired_tumornormal_calling(install_test_files, data_dir, global_config):
    # using MuTect, VarScan, FreeBayes
    run_config = os.path.join(data_dir, 'run_info-cancer.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, run_config])


@pytest.mark.cancer
@pytest.mark.cancerpanel
@pytest.mark.install_required
def test_cancer_calling_with_and_without_normal(install_test_files, data_dir, global_config):
    # requires MuTect and GATK
    run_config = os.path.join(data_dir, 'run_info-cancer2.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, run_config])


@pytest.mark.cancer
@pytest.mark.cancerprecall
@pytest.mark.install_required
def test_cancer_with_precalled_inputs(install_test_files, data_dir, global_config):
    run_config = os.path.join(data_dir, 'run_info-cancer3.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, run_config])


@pytest.mark.speed1
@pytest.mark.template
def test_create_project_template(install_test_files, data_dir, work_dir):
    fc_dir = os.path.join(data_dir, os.pardir, '100326_FC6107FAAXX')
    subprocess.check_call(['bcbio_nextgen.py', '-w', 'template', '--only-metadata',
                           'freebayes-variant', os.path.join(fc_dir, '100326.csv'),
                           os.path.join(fc_dir, '7_100326_FC6107FAAXX_1_fastq.txt'),
                           os.path.join(fc_dir, '7_100326_FC6107FAAXX_2_fastq.txt'),
                           os.path.join(fc_dir, '8_100326_FC6107FAAXX.bam')])


@pytest.mark.joint
@pytest.mark.install_required
def test_joint_calling(install_test_files, data_dir, global_config):
    run_config = os.path.join(data_dir, 'run_info-joint.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, run_config])


@pytest.mark.umibarcode
@pytest.mark.install_required
def test_umi_with_bam_inputs(install_test_files, data_dir, global_config):
    run_config = os.path.join(data_dir, 'run_info-umi.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, run_config])


@pytest.mark.hla
def test_hla_typing_with_optitype(install_test_files, data_dir, work_dir):
    from bcbio.hla import optitype
    hla_dir = os.path.join(data_dir, os.pardir, '100326_FC6107FAAXX', 'hla')
    data = {"dirs": {"work": work_dir},
            "rgnames": {"sample": "test"},
            "config": {},
            "hla": {"fastq": glob.glob(os.path.join(hla_dir, "*"))}}
    out = optitype.run(data)
    with open(out["hla"]["call_file"]) as in_handle:
        header = in_handle.readline().strip().split(",")
        hla_a = dict(zip(header, in_handle.readline().strip().split(",")))
        assert hla_a["alleles"] == "HLA-A*11:01;HLA-A*24:02", hla_a
