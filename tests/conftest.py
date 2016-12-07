import collections
import contextlib
import os
import shutil
import subprocess

import pytest
import yaml

from bcbio.pipeline.config_utils import load_system_config

OUTPUT_DIR = "test_automated_output"


def default_workdir():
    return os.path.join(os.path.dirname(__file__), OUTPUT_DIR)


@pytest.fixture
def data_dir():
    return os.path.join(os.path.dirname(__file__), "data", "automated")


@contextlib.contextmanager
def make_workdir():
    remove_old_dir = True
    # Specify workdir though env var, in case tests have to run not in the
    # default location (e.g. to run tests on a  mounted FS)
    dirname = os.environ.get('BCBIO_WORKDIR', default_workdir())
    if remove_old_dir:
        if os.path.exists(dirname):
            shutil.rmtree(dirname)
        os.makedirs(dirname)
    orig_dir = os.getcwd()
    try:
        os.chdir(dirname)
        yield dirname
    finally:
        os.chdir(orig_dir)


@pytest.yield_fixture
def workdir():
    with make_workdir() as wd:
        yield wd


def get_post_process_yaml(data_dir, workdir):
    """Prepare a bcbio_system YAML file pointing to test data.
    """
    try:
        from bcbiovm.docker.defaults import get_datadir
        datadir = data_dir or get_datadir()
        sys_conf_file = os.path.join(datadir, "galaxy", "bcbio_system.yaml")
        system = sys_conf_file if datadir else None
    except ImportError:
        system = None
    if system is None or not os.path.exists(system):
        try:
            _, system = load_system_config(
                config_file="bcbio_system.yaml", work_dir=workdir)
        except ValueError:
            system = None
    if system is None or not os.path.exists(system):
        system = os.path.join(data_dir, "post_process-sample.yaml")
    # create local config pointing to reduced genomes
    test_system = os.path.join(workdir, "bcbio_system.yaml")
    with open(system) as in_handle:
        config = yaml.load(in_handle)
        config["galaxy_config"] = os.path.join(data_dir, "universe_wsgi.ini")
        with open(test_system, "w") as out_handle:
            yaml.dump(config, out_handle)
    return test_system


@pytest.fixture
def install_test_files(data_dir):
    """Download required sequence and reference files.
    """
    DlInfo = collections.namedtuple("DlInfo", "fname dirname version")
    download_data = [
        DlInfo("110106_FC70BUKAAXX.tar.gz", None, None),
        DlInfo("genomes_automated_test.tar.gz", "genomes", 31),
        DlInfo("110907_ERP000591.tar.gz", None, None),
        DlInfo("100326_FC6107FAAXX.tar.gz", None, 11),
        DlInfo("tcga_benchmark.tar.gz", None, 3),
        DlInfo("singlecell-rnaseq-test-data.tar.gz", "Harvard-inDrop", 1)
    ]
    for dl in download_data:
        url = "http://chapmanb.s3.amazonaws.com/{fname}".format(fname=dl.fname)
        dirname = os.path.join(
            data_dir, os.pardir,
            dl.fname.replace(".tar.gz", "") if dl.dirname is None
            else dl.dirname
        )
        if os.path.exists(dirname) and dl.version is not None:
            version_file = os.path.join(dirname, "VERSION")
            is_old = True
            if os.path.exists(version_file):
                with open(version_file) as in_handle:
                    version = int(in_handle.read())
                is_old = version < dl.version
            if is_old:
                shutil.rmtree(dirname)
        if not os.path.exists(dirname):
            _download_to_dir(url, dirname)


def _download_to_dir(url, dirname):
        print(dirname)
        cl = ["wget", url]
        subprocess.check_call(cl)
        cl = ["tar", "-xzvpf", os.path.basename(url)]
        subprocess.check_call(cl)
        shutil.move(os.path.basename(dirname), dirname)
        os.remove(os.path.basename(url))
