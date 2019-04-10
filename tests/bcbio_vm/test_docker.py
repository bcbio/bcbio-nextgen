import os
import shutil
import subprocess

import pytest

from bcbio import utils

from tests.conftest import (make_workdir, get_post_process_yaml, install_cwl_test_files)

@pytest.mark.docker
@pytest.mark.docker_multicore
def test_docker(install_test_files, data_dir):
    """Run an analysis with code and tools inside a docker container.

    Requires https://github.com/bcbio/bcbio-nextgen-vm
    """
    with make_workdir() as workdir:
        cl = [
            "bcbio_vm.py",
            "--datadir=%s" % data_dir,
            "run",
            "--image=quay.io/bcbio/bcbio-vc",
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

    Requires https://github.com/bcbio/bcbio-nextgen-vm
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

    Requires https://github.com/bcbio/bcbio-nextgen-vm
    """
    @pytest.mark.cwl
    @pytest.mark.cwl_docker
    @pytest.mark.cwl_docker_somatic
    def test_2_cwl_docker_somatic(self, data_dir):
        """CWL: run a somatic workflow using Docker.
        """
        with install_cwl_test_files(data_dir) as workdir:
            with utils.chdir(os.path.join(workdir, "somatic")):
                cl = ["bash", "./run_generate_cwl.sh"]
                subprocess.check_call(cl)
                if os.path.exists("cromwell_work"):
                    shutil.rmtree("cromwell_work")
                cl = ["bcbio_vm.py", "cwlrun", "cromwell", "somatic-workflow"]
                subprocess.check_call(cl)

    @pytest.mark.cwl
    @pytest.mark.cwl_docker
    @pytest.mark.cwl_docker_joint
    def test_3_cwl_docker_joint(self, data_dir):
        """CWL: run a gVCF based joint-calling workflow using Docker.
        """
        with install_cwl_test_files(data_dir) as workdir:
            with utils.chdir(os.path.join(workdir, "gvcf_joint")):
                cl = ["bash", "./run_generate_cwl.sh"]
                subprocess.check_call(cl)
                if os.path.exists("cromwell_work"):
                    shutil.rmtree("cromwell_work")
                cl = ["bash", "./run_cromwell.sh"]
                subprocess.check_call(cl)

    @pytest.mark.cwl
    @pytest.mark.cwl_docker_rnaseq
    def test_4_cwl_rnaseq(self, install_test_files, data_dir):
        """CWL: run a RNA-seq workflow using Docker.
        """
        with install_cwl_test_files(data_dir) as workdir:
            with utils.chdir(os.path.join(workdir, "rnaseq")):
                cl = ["bash", "./run_generate_cwl.sh"]
                if os.path.exists("cromwell_work"):
                    shutil.rmtree("cromwell_work")
                cl = ["bcbio_vm.py", "cwlrun", "cromwell", "rnaseq-workflow"]
                subprocess.check_call(cl)

    @pytest.mark.cwl
    @pytest.mark.cwl_local
    @pytest.mark.install_required
    def test_1_cwl_local(self, install_test_files, data_dir):
        """CWL: prepare somatic workflow and run on local installation.
        """
        with install_cwl_test_files(data_dir) as workdir:
            with utils.chdir(os.path.join(workdir, "somatic")):
                cl = ["bash", "./run_generate_cwl.sh"]
                subprocess.check_call(cl)
                if os.path.exists("cwltool_work"):
                    shutil.rmtree("cwltool_work")
                cl = ["bash", "./run_cwltool.sh"]
                subprocess.check_call(cl)

    @pytest.mark.cwl
    @pytest.mark.cwl_arvados
    def test_5_cwl_arvados(self, install_test_files, data_dir):
        """CWL: prepare workflow for running on Arvados.

        Requires ARVADOS_API_HOST and ARVADOS_API_TOKEN set
        """
        if os.environ.get("ARVADOS_API_HOST") and os.environ.get("ARVADOS_API_TOKEN"):
            with install_cwl_test_files(data_dir) as workdir:
                with utils.chdir(os.path.join(workdir, "arvados")):
                    cl = ["bash", "./run_generate_cwl.sh"]
                    subprocess.check_call(cl)
