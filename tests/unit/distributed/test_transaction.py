import pytest
import mock

from bcbio.distributed import transaction
from bcbio.distributed.transaction import tx_tmpdir
from bcbio.distributed.transaction import file_transaction
from bcbio.distributed.transaction import _get_base_tmpdir
from bcbio.distributed.transaction import _flatten_plus_safe
from bcbio.distributed.transaction import _move_file_with_sizecheck
from tests.unit.conftest import DummyCM, DummyTxTmpdir


CWD = 'TEST_CWD'
CONFIG = {'a': 1}
TMP = '/tmp'
TMPDIR = 'TEST_TMPDIR'


class DummyFlattenPlusSafe(DummyCM):
    value = (['foo'], ['bar'])

    def __iter__(self):
        for v in self.value:
            yield v


@pytest.yield_fixture
def mock_flatten(mocker):
    yield mocker.patch(
        'bcbio.distributed.transaction._flatten_plus_safe',
        side_effect=DummyFlattenPlusSafe
    )


@pytest.yield_fixture
def mock_io(mocker):
    mocker.patch('bcbio.distributed.transaction.open')
    mocker.patch('bcbio.distributed.transaction.os.path.isdir')
    mocker.patch('bcbio.distributed.transaction.os.path.isfile')
    mocker.patch('bcbio.distributed.transaction.os.path.exists')
    mocker.patch('bcbio.distributed.transaction.shutil')
    mocker.patch(
        'bcbio.distributed.transaction.tempfile.mkdtemp',
        return_value=TMPDIR)
    mocker.patch('bcbio.distributed.transaction.utils')
    mocker.patch(
        'bcbio.distributed.transaction.os.getcwd',
        return_value=CWD
    )
    yield None


class TestTxTmpdir(object):

    def test_gets_base_tmpdir_name_from_config_or_cwd(self, mock_io, mocker):
        mocker.patch('bcbio.distributed.transaction._get_base_tmpdir')
        data  = mock.Mock()
        with tx_tmpdir(data):
            pass
        cwd = transaction.os.getcwd.return_value
        transaction._get_base_tmpdir.assert_called_once_with(
            data, cwd)
        base_tmpdir = transaction._get_base_tmpdir.return_value
        transaction.utils.get_abspath.assert_called_once_with(base_tmpdir)

    def test_makes_base_tmp_dir(self, mock_io):
        """"
        Test that tx_tmpdir creates a base temporary directory
        """
        with tx_tmpdir(None):
            pass
        transaction.utils.safe_makedir.assert_called_once_with(
            transaction.utils.get_abspath.return_value)

    def test_makes_unique_tmp_dir(self, mock_io):
        """Test that tx_tmpdir creates a tmp dir unique name
        using `tempfile.mkdtemp` inside the base dir."""
        with tx_tmpdir(None):
            pass
        transaction.tempfile.mkdtemp.assert_called_once_with(
            dir=transaction.utils.get_abspath.return_value)

    def test_yields_tmp_dir(self, mock_io):
        """Test that tx_tmpdir yields a path to the created directory."""
        expected = transaction.tempfile.mkdtemp.return_value
        with tx_tmpdir() as tmp_dir:
            assert tmp_dir == expected

    def test_rmtree_not_called_if_remove_is_false(self, mock_io):
        with tx_tmpdir(remove=False):
            pass
        assert not transaction.utils.remove_safe.called

    def test_rmtree_called_if_remove_is_true(self, mock_io):
        transaction.tempfile.mkdtemp.return_value = 'foo'
        with tx_tmpdir(remove=True):
            pass
        transaction.utils.remove_safe.assert_called_once_with('foo')

    def test_create_tmpdir_in_a_specified_base_dir(self, mock_io):
        with tx_tmpdir(base_dir='somedir'):
            pass
        transaction.utils.get_abspath.assert_called_once_with(
            'somedir/bcbiotx')
        transaction.utils.safe_makedir.assert_called_once_with(
            transaction.utils.get_abspath.return_value)


class TestGetConfigTmpdir(object):
    """Test that base tmpdir is extracted from config properly"""
    def test_get_config_tmpdir__from_config(self):
        TMPDIR = 'TEST_TMP_DIR'
        config = {
            'config': {
                'resources': {
                    'tmp': {'dir': TMPDIR}
                }
            }
        }
        expected = TMPDIR
        result = _get_base_tmpdir(config, CWD)
        assert result == expected

    def test_get_config_tmpdir__from_resources(self):
        config = {
            'resources': {
                'tmp': {'dir': 'TEST_TMP_DIR'}
            }
        }
        expected = 'TEST_TMP_DIR'
        result = _get_base_tmpdir(config, CWD)
        assert result == expected

    def test_get_config_tmpdir__no_data(self):
        result = _get_base_tmpdir(None, CWD)
        assert result == '%s/bcbiotx' % CWD


class TestFlattenPlusSafe(object):
    """Tests the logic of handling arguments passed to file_transaction
    in different forms.
    """
    @pytest.yield_fixture(autouse=True)
    def mock_tx_tmpdir(self, mocker):
        yield mocker.patch(
            'bcbio.distributed.transaction.tx_tmpdir',
            side_effect=DummyTxTmpdir
        )

    @pytest.mark.parametrize(('args', 'exp_tx_args'), [
        (('/path/to/somefile',), (None,)),
        ((CONFIG, '/path/to/somefile'), (CONFIG,)),
        ((CONFIG, ['/path/to/somefile']), (CONFIG,)),
        ((CONFIG, '/path/to/somefile', '/otherpath/to/file'), (CONFIG,))
    ])
    def test_calls_tx_tmpdir(self, args, exp_tx_args):
        with _flatten_plus_safe(args) as (result_tx, result_safe):
            pass
        transaction.tx_tmpdir.assert_called_once_with(*exp_tx_args)

    @pytest.mark.parametrize(('args', 'expected_tx'), [
        (('/path/to/somefile',), ['foo/somefile']),
        ((CONFIG, '/path/to/somefile'), ['foo/somefile']),
        ((CONFIG, ['/path/to/somefile']), ['foo/somefile']),
        (
            (CONFIG, '/path/to/somefile', '/otherpath/to/otherfile'),
            ['foo/somefile', 'foo/otherfile'],
        )]
    )
    def test_creates_path_to_tx_file_in_tmp_dir(self, args, expected_tx):
        with _flatten_plus_safe(args) as (result_tx, _):
            assert result_tx == expected_tx

    @pytest.mark.parametrize(('args', 'expected_safe'), [
        (('/path/to/somefile',), ['/path/to/somefile']),
        ((CONFIG, '/path/to/somefile'), ['/path/to/somefile']),
        ((CONFIG, ['/path/to/somefile']), ['/path/to/somefile']),
        (
            (CONFIG, '/path/to/somefile', '/otherpath/to/otherfile'),
            ['/path/to/somefile', '/otherpath/to/otherfile'],
        )]
    )
    def test_returns_original_filesnames(self, args, expected_safe):
        with _flatten_plus_safe(args) as (_, result_safe):
            assert result_safe == expected_safe


class TestMoveWithSizeCheck(object):
    def test_moves_files(self, mock_io):
        _move_file_with_sizecheck('foo', 'bar')
        transaction.shutil.move.assert_called_once_with('foo', 'bar')

    def test_fails_if_sizes_arent_equal(self, mock_io):
        transaction.utils.get_size.side_effect = lambda x: x
        with pytest.raises(AssertionError):
            _move_file_with_sizecheck('foo', 'bar')
        assert not transaction.utils.remove_safe.called

    def test_creates_flag_file(self, mock_io):
        expected_flag_file = 'bar.bcbiotmp'
        _move_file_with_sizecheck('foo', 'bar')
        transaction.open.assert_called_once_with(expected_flag_file, 'wb')
        transaction.utils.remove_safe.assert_called_once_with(
            expected_flag_file)


class TestFileTransaction(object):

    def test_yields_path_to_tmpfile(self, mock_io):
        original = '/some/path'
        expected_tmp = '%s/path' % TMPDIR
        with file_transaction(CONFIG, original) as tmp_path:
            assert tmp_path == expected_tmp

    def test_yields_tuple_of_paths_to_tmpfiles(self, mock_io):
        original = ('/some/path', '/some/otherpath')
        expected_tmp = ('%s/path' % TMPDIR, '%s/otherpath' % TMPDIR)
        with file_transaction(CONFIG, *original) as tmp_path:
            assert tmp_path == expected_tmp

    def test_moves_tmp_file_to_original_locations(self, mock_io):
        original = '/some/path'
        transaction.os.path.exists.return_value = True
        with file_transaction(CONFIG, original) as tmp_path:
            pass
        transaction.shutil.move.assert_called_once_with(tmp_path, original)

    def test_removes_tmpdir(self, mock_io):
        with file_transaction(CONFIG, '/some/path'):
            pass
        transaction.utils.remove_safe.assert_called_with(TMPDIR)

    def test_doesnt_try_to_move_tmp_file_if_it_doesnt_exist(self, mock_io):
        transaction.os.path.exists.return_value = False
        with file_transaction(CONFIG, '/some/path'):
            pass
        assert not transaction.shutil.move.called
