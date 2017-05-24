import os
import subprocess

import pytest
from tests.conftest import make_workdir
from tests.conftest import get_post_process_yaml


@pytest.mark.S3
@pytest.mark.install_required
def test_fusion(install_test_files, data_dir):
    """Run an RNA-seq analysis and test fusion genes, with human-mouse
    disambiguation.
    Requires minikraken database.
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, os.pardir, "test_fusion"),
              os.path.join(data_dir, "run_info-fusion_S3.yaml")]
        subprocess.check_call(cl)


@pytest.mark.S3
@pytest.mark.install_required
def test_variantcall_1(install_test_files, data_dir):
    """Test variant calling with disambiguation.
    Requires minikraken database.
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, os.pardir, "100326_FC6107FAAXX"),
              os.path.join(data_dir, "run_info-variantcall_S3_1.yaml")]
        subprocess.check_call(cl)


@pytest.mark.S3
@pytest.mark.install_required
def test_variantcall_2(install_test_files, data_dir):
    """Test variant calling with disambiguation.
    Requires minikraken database.
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, os.pardir, "100326_FC6107FAAXX"),
              os.path.join(data_dir, "run_info-variantcall_S3_2.yaml")]
        subprocess.check_call(cl)
