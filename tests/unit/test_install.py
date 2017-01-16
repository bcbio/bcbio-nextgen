import mock
import pytest

import subprocess
from bcbio.install import CondaAPI


class TestCondaAPI(object):

    @pytest.yield_fixture
    def sp(self, mocker):
        yield mocker.patch('bcbio.install.subprocess.check_output')

    @pytest.fixture
    def api(self):
        return CondaAPI()

    def test_create_env_returns_prefix_to_env(self, api, sp):
        sp.return_value = '{"actions": {"PREFIX": "/path/to/env"}}'
        result = api.create_env('test')
        assert result == '/path/to/env'

    def test_create_env_calls_subprocess_check_output(self, api, sp):
        sp.return_value = '{"actions": {"PREFIX": "/path/to/env"}}'
        name = 'test_name'
        api.create_env(name)
        expected = [
            'conda', 'create',
            '--name', name,
            '-c', 'bioconda', '-c', 'r', '-c', 'conda-forge',
            '--json',
        ]
        sp.assert_called_once_with(expected)

    def test_create_env_parses_error_message_if_already_exists(self, api, sp):
        msg = '{"message": "Value error: prefix already exists: /path/to/env"}'
        sp.side_effect = subprocess.CalledProcessError(
            output=msg, returncode=1, cmd=mock.Mock())
        result = api.create_env('test_name')
        assert result == '/path/to/env'

    def test_raises_subprocess_error_if_it_is_unknown(self, api, sp):
        msg = '{"message": "bla bla bla"}'
        sp.side_effect = subprocess.CalledProcessError(
            output=msg, returncode=1, cmd=mock.Mock())
        with pytest.raises(subprocess.CalledProcessError):
            api.create_env('test_name')

    def test_get_latest_version_calls_subprocess_check_output(self, sp, api):
        pkg = 'test_pkg'
        sp.return_value = '{"%s": [{"version": "1.0.0"}]}' % pkg
        api.get_latest_version(pkg)
        sp.assert_called_once_with([
            'conda', 'search', pkg,
            '-c', 'bioconda', '-c', 'r', '-c', 'conda-forge',
            '--json',
        ])

    def test_get_lates_version_returns_package_version(self, api, sp):
        pkg = 'test_pkg'
        sp.return_value = \
            '{"%s": [{"version": "1.0"},{"version": "1.1"}]}' % pkg
        result = api.get_latest_version(pkg)
        assert result == "1.1"

    def test_install_package_latest_version(self, api, sp):
        pkg = 'foo'
        sp.return_value = '{"success": true}'
        api.install_package(pkg)
        sp.assert_called_once_with([
            'conda', 'install',
            '--name', 'root',
            pkg,
            '--quiet',
            '-c', 'bioconda', '-c', 'r', '-c', 'conda-forge',
            '--json',
        ])

    def test_install_package_specific_version(self, api, sp):
        pkg = 'foo'
        version = '1.2.3'
        sp.return_value = '{"success": true}'
        api.install_package(pkg, version=version)
        sp.assert_called_once_with([
            'conda', 'install',
            '--name', 'root',
            '%s=%s' % (pkg, version),
            '--quiet',
            '-c', 'bioconda', '-c', 'r', '-c', 'conda-forge',
            '--json',
        ])

    def test_install_package_into_specific_environment(self, api, sp):
        pkg = 'foo'
        env_name = 'test_env'
        sp.return_value = '{"success": true}'
        api.install_package(pkg, env_name=env_name)
        sp.assert_called_once_with([
            'conda', 'install',
            '--name', env_name,
            pkg,
            '--quiet',
            '-c', 'bioconda', '-c', 'r', '-c', 'conda-forge',
            '--json',
        ])

    def test_install_package_raises_if_unsuccessful(self, api, sp):
        sp.return_value = '{"success": false}'
        with pytest.raises(RuntimeError):
            api.install_package('test_pkg')
