import os
import subprocess

from tests.conftest import make_workdir
from tests.conftest import get_post_process_yaml


def test_fusion(install_test_files, data_dir):
    """Run an RNA-seq analysis and test fusion genes, with human-mouse
    disambiguation.
    """
    with make_workdir() as workdir:
        cl = ["bcbio_nextgen.py",
              get_post_process_yaml(data_dir, workdir),
              os.path.join(data_dir, os.pardir, "test_fusion"),
              os.path.join(data_dir, "run_info-fusion_S3.yaml")]
        subprocess.check_call(cl)
