import os
import subprocess

import pytest


@pytest.mark.S3
@pytest.mark.install_required
@pytest.mark.xfail(reason='https://github.com/bcbio/bcbio-nextgen/issues/3216', run=False)
def test_fusion(install_test_files, data_dir, global_config):
    """Run an RNA-seq analysis and test fusion genes, with human-mouse disambiguation.
    Requires minikraken database.
    """
    fc_dir = os.path.join(data_dir, os.pardir, 'test_fusion')
    run_config = os.path.join(data_dir, 'run_info-fusion_S3.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])


@pytest.mark.S3
@pytest.mark.install_required
@pytest.mark.xfail(reason='https://github.com/bcbio/bcbio-nextgen/issues/3217', run=False)
def test_variantcall_1(install_test_files, data_dir, global_config):
    """Test variant calling with disambiguation.
    Requires minikraken database.
    """
    fc_dir = os.path.join(data_dir, os.pardir, '100326_FC6107FAAXX')
    run_config = os.path.join(data_dir, 'run_info-variantcall_S3_1.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])


@pytest.mark.S3
@pytest.mark.install_required
@pytest.mark.xfail(reason='https://github.com/bcbio/bcbio-nextgen/issues/3218', run=False)
def test_variantcall_2(install_test_files, data_dir, global_config):
    """Test variant calling with disambiguation.
    Requires minikraken database.
    """
    fc_dir = os.path.join(data_dir, os.pardir, '100326_FC6107FAAXX')
    run_config = os.path.join(data_dir, 'run_info-variantcall_S3_2.yaml')
    subprocess.check_call(['bcbio_nextgen.py', global_config, fc_dir, run_config])
