"""Pytest fixtures and test helper functions"""

import collections
import contextlib
from datetime import datetime
import io
import os
import shutil
import subprocess
import tarfile
import tempfile

import pytest
import requests
import yaml

from bcbio.pipeline.config_utils import load_system_config

if os.environ.get("BCBIO_TEST_DIR"):
    BCBIO_TEST_DIR = os.environ.get("BCBIO_TEST_DIR")
else:
    BCBIO_TEST_DIR = tempfile.TemporaryDirectory(prefix="bcbio_").name  # /tmp/bcbio

def pytest_addoption(parser):
    parser.addoption('--keep-test-dir', action='store_true', default=False,
                     help='Preserve test output directory after each test')


@pytest.fixture(scope='session')
def test_dir(pytestconfig):
    os.makedirs(BCBIO_TEST_DIR, exist_ok=True)
    yield
    if not pytestconfig.getoption('--keep-test-dir'):
        shutil.rmtree(BCBIO_TEST_DIR)


@pytest.fixture(scope='session')
def data_dir(test_dir):
    # workaround for hardcoded data file paths in test run config files
    test_data_dir = os.path.join(BCBIO_TEST_DIR, 'data')  # /tmp/bcbio/data
    with contextlib.suppress(FileExistsError):
        os.symlink(os.path.join(os.path.dirname(__file__), 'data'), test_data_dir)
    return os.path.join(test_data_dir, 'automated')  # /tmp/bcbio/data/automated


@pytest.fixture
def work_dir(pytestconfig):
    """Provide and manage output directory for tests"""
    test_output_dir = os.path.join(BCBIO_TEST_DIR, 'test_automated_output')
    os.makedirs(test_output_dir, exist_ok=True)
    original_dir = os.getcwd()
    os.chdir(test_output_dir)
    yield test_output_dir
    os.chdir(original_dir)
    if not pytestconfig.getoption('--keep-test-dir'):
        shutil.rmtree(test_output_dir)


@pytest.fixture
def global_config(data_dir, work_dir):
    """Prepare a bcbio_system YAML file pointing to test data"""
    system = _get_bcbio_system(work_dir, data_dir)
    # create local config pointing to reduced genomes
    test_system = os.path.join(work_dir, 'bcbio_system.yaml')
    with open(system) as in_handle:
        config = yaml.safe_load(in_handle)
        config["galaxy_config"] = os.path.join(data_dir, "universe_wsgi.ini")
        with open(test_system, "w") as out_handle:
            yaml.dump(config, out_handle)
    return test_system


@contextlib.contextmanager
def install_cwl_test_files():
    orig_dir = os.getcwd()
    url = "https://github.com/bcbio/test_bcbio_cwl/archive/master.tar.gz"
    dirname = os.path.join(BCBIO_TEST_DIR, 'test_bcbio_cwl-master')
    if os.path.exists(dirname):
        # check for updated commits if the directory exists
        ctime = os.path.getctime(os.path.join(dirname, "README.md"))
        dtime = datetime.fromtimestamp(ctime).isoformat()
        r = requests.get("https://api.github.com/repos/bcbio/test_bcbio_cwl/commits?since=%s" % dtime).json()
        if len(r) > 0:
            shutil.rmtree(dirname)
    try:
        if not os.path.exists(dirname):
            print("Downloading CWL test directory: %s" % url)
            os.makedirs(dirname)
            os.chdir(os.path.dirname(dirname))
            r = requests.get(url)
            tf = tarfile.open(fileobj=io.BytesIO(r.content), mode='r|gz')
            tf.extractall()
        os.chdir(dirname)
        yield dirname
    finally:
        os.chdir(orig_dir)


def _get_bcbio_system(workdir, data_dir):
    system = _get_bcbiovm_config(data_dir)
    if _config_is_invalid(system):
        system = _get_system_config(workdir, system)
    if _config_is_invalid(system):
        system = os.path.join(data_dir, "post_process-sample.yaml")
    return system


def _get_bcbiovm_config(data_dir):
    try:
        from bcbiovm.docker.defaults import get_datadir
        datadir = data_dir or get_datadir()
        sys_conf_file = os.path.join(datadir, "galaxy", "bcbio_system.yaml")
        system = sys_conf_file if datadir else None
    except ImportError:
        system = None
    return system


def _config_is_invalid(config_fname):
    return config_fname is None or not os.path.exists(config_fname)


def _get_system_config(work_dir, system):
    try:
        _, system = load_system_config(
            config_file="bcbio_system.yaml",
            work_dir=work_dir
        )
    except ValueError:
        system = None
    return system


@pytest.fixture
def install_test_files(data_dir):
    """Download required sequence and reference files."""
    DlInfo = collections.namedtuple("DlInfo", "fname dirname version")
    download_data = [
        DlInfo("110106_FC70BUKAAXX.tar.gz", None, None),
        DlInfo("genomes_automated_test.tar.gz", "genomes", 34),
        DlInfo("110907_ERP000591.tar.gz", None, None),
        DlInfo("100326_FC6107FAAXX.tar.gz", None, 12),
        DlInfo("tcga_benchmark.tar.gz", None, 3),
        DlInfo("singlecell-rnaseq-test-data.tar.gz", "Harvard-inDrop", 1)
    ]
    for dl in download_data:
        url = f"https://bcbio-nextgen.s3.amazonaws.com/test_data/{dl}"
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
    subprocess.check_call(['wget', '--progress=dot:giga', url])
    subprocess.check_call(['tar', '-xzvpf', os.path.basename(url)])
    shutil.move(os.path.basename(dirname), dirname)
    os.remove(os.path.basename(url))
