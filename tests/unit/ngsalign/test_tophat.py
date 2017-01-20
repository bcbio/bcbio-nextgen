from bcbio.ngsalign import tophat

import mock
import pytest


@pytest.yield_fixture
def mock_should_run(mocker):
    yield mocker.patch(
        'bcbio.ngsalign.tophat.config_utils.should_run_fusion')


def test_should_run_fusion(mock_should_run):
    config = mock.MagicMock()
    result = tophat._should_run_fusion(config)
    mock_should_run.assert_called_once_with('tophat', config)
    assert result == mock_should_run.return_value


def test_set_fusion_mode_when_should_run_fusion(mock_should_run):
    mock_should_run.return_value = True
    opts = {}
    tophat._set_fusion_mode(opts, mock.Mock())
    assert opts['fusion-search'] is True


def test_set_fusion_mode_when_should_not_run_fusion(mock_should_run):
    mock_should_run.return_value = False
    opts = {}
    tophat._set_fusion_mode(opts, mock.Mock())
    assert opts == {}
