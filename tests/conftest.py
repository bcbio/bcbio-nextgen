import collections
import contextlib
from datetime import datetime
import io
import os
import shutil
import subprocess
import tarfile

import pytest
import requests
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
    system = _get_bcbio_system(workdir, data_dir)
    # create local config pointing to reduced genomes
    test_system = os.path.join(workdir, "bcbio_system.yaml")
    with open(system) as in_handle:
        config = yaml.safe_load(in_handle)
        config["galaxy_config"] = os.path.join(data_dir, "universe_wsgi.ini")
        with open(test_system, "w") as out_handle:
            yaml.dump(config, out_handle)
    return test_system

@contextlib.contextmanager
def install_cwl_test_files(data_dir):
    orig_dir = os.getcwd()
    url = "https://github.com/bcbio/test_bcbio_cwl/archive/master.tar.gz"
    dirname = os.path.normpath(os.path.join(data_dir, os.pardir, "test_bcbio_cwl-master"))
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
    """Download required sequence and reference files.
    """
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
        cl = ["wget", url]
        subprocess.check_call(cl)
        cl = ["tar", "-xzvpf", os.path.basename(url)]
        subprocess.check_call(cl)
        shutil.move(os.path.basename(dirname), dirname)
        os.remove(os.path.basename(url))
