import mock
import pytest

from bcbio import utils


@pytest.yield_fixture
def bcbio_env(mocker):
    env = mocker.patch('bcbio.utils.get_bcbio_env')
    env.return_value = {'PATH': '/usr/bin:/bin'}
    yield env


def test_get_ericscript_env_if_executable_not_found(bcbio_env):
    test_conda_prefix = '/path/to/conda/env'
    env = utils.get_ericscript_env(test_conda_prefix)
    assert env['PATH'].startswith(test_conda_prefix)


def test_get_ericscript_env_raises_runtime_error_if_no_env_prefix(bcbio_env):
    with pytest.raises(RuntimeError):
        utils.get_ericscript_env(None)
