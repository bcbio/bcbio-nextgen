import os
import subprocess

import pytest

from tests.conftest import make_workdir
from tests.conftest import get_post_process_yaml

@pytest.mark.docker
def test_docker(install_test_files, data_dir):
    """Run an analysis with code and tools inside a docker container.

    Requires https://github.com/chapmanb/bcbio-nextgen-vm
    """
    with make_workdir() as workdir:
        cl = [
            "bcbio_vm.py",
            "--datadir=%s" % data_dir,
            "run",
            "--systemconfig=%s" % get_post_process_yaml(data_dir, workdir),
            "--fcdir=%s" % os.path.join(
                data_dir, os.pardir, "100326_FC6107FAAXX"),
            os.path.join(data_dir, "run_info-bam.yaml")
        ]
        subprocess.check_call(cl)

@pytest.mark.docker
@pytest.mark.docker_ipython
def test_docker_ipython(install_test_files, data_dir):
    """Run an analysis with code and tools inside a docker container,
    driven via IPython.

    Requires https://github.com/chapmanb/bcbio-nextgen-vm
    """
    with make_workdir() as workdir:
        cl = [
            "bcbio_vm.py",
            "--datadir=%s" % data_dir,
            "ipython",
            "--systemconfig=%s" % get_post_process_yaml(data_dir, workdir),
            "--fcdir=%s" % os.path.join(
                data_dir, os.pardir, "100326_FC6107FAAXX"),
            os.path.join(data_dir, "run_info-bam.yaml"),
            "lsf", "localrun"
        ]
        subprocess.check_call(cl)


class TestCWL():
    """ Run simple CWL workflows.

    Requires https://github.com/chapmanb/bcbio-nextgen-vm
    """
    @pytest.mark.cwl_docker
    @pytest.mark.cwl
    def test_2_cwl_docker(install_test_files, data_dir):
        """Create a common workflow language description and run on a
        Docker installation.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_vm.py", "cwl", "../data/automated/run_info-cwl.yaml",
                  "--systemconfig", get_post_process_yaml(data_dir, workdir)]
            subprocess.check_call(cl)
            cl = ["bcbio_vm.py", "cwlrun", "cwltool", "run_info-cwl-workflow"]
            subprocess.check_call(cl)
            print
            print "To run with a CWL tool, cd test_automated_output and:"
            print " ".join(cl)

    @pytest.mark.speed2
    @pytest.mark.cwl
    @pytest.mark.cwl_local
    @pytest.mark.install_required
    def test_1_cwl_local(self, install_test_files, data_dir):
        """Create a common workflow language description and run on local installation.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_vm.py", "cwl", "../data/automated/run_info-cwl.yaml",
                  "--systemconfig", get_post_process_yaml(data_dir, workdir)]
            subprocess.check_call(cl)
            cl = ["bcbio_vm.py", "cwlrun", "cwltool", "run_info-cwl-workflow",
                  "--no-container"]
            subprocess.check_call(cl)
            print
            print "To run with a CWL tool, cd test_automated_output and:"
            print " ".join(cl)
