import mock
import pytest

from bcbio import utils


@pytest.yield_fixture
def environ(mocker):
    environ = mocker.patch('bcbio.utils.os.environ.copy')
    environ.return_value = {'PATH': '/usr/bin:/bin'}


def test_get_ericscript_env(environ):
    env = utils.get_ericscript_env(mock.Mock())
    expected = '/usr/local/share/bcbio-nextgen/anaconda/envs/ericscript/bin'
    assert env['PATH'] == '%s:/usr/bin:/bin' % expected
