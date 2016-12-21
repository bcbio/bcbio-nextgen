import mock
import pytest

from bcbio import utils

TEST_CONDA_ENV = '/path/to/conda/env'


@pytest.yield_fixture
def bcbio_env(mocker):
    env = mocker.patch('bcbio.utils.get_bcbio_env')
    env.return_value = {'PATH': '/usr/bin:/bin'}
    yield env


@pytest.yield_fixture
def which(mocker):
    yield mock.patch('bcbio.utils.which')


@pytest.fixture
def config():
    return {'config': {'resources': {'ericscript': {'env': TEST_CONDA_ENV}}}}


def test_get_ericscript_env_if_executable_not_found(which, bcbio_env, config):
    which.return_value = None
    env = utils.get_ericscript_env(config)
    assert env['PATH'].startswith(TEST_CONDA_ENV)


def test_ericscript_env_if_exec_found_in_bcbio_env(which, bcbio_env, config):
    which.return_value = '/some/path'
    env = utils.get_ericscript_env(config)
    assert env == bcbio_env.return_value


def test_get_ericscript_env_raises_runtime_error_when_install_not_found(
        which, bcbio_env):
    which.return_value = None
    config = {'config': {'resources': {}}}
    with pytest.raises(RuntimeError):
        utils.get_ericscript_env(config)
