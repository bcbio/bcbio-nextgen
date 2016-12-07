from bcbio.ngsalign import tophat

import mock
import pytest


@pytest.yield_fixture
def mock_should_run(mocker):
    yield mocker.patch('bcbio.ngsalign.tophat._should_run_fusion')


@pytest.mark.parametrize('config', [
    {'algorithm': {'fusion_mode': False}},
    {'config': {'algorithm': {'fusion_mode': False}}},
])
def test_should_not_run_fusion_when_fusion_mode_false(config):
    result = tophat._should_run_fusion(config)
    assert result is False


@pytest.mark.parametrize('config', [
    {'algorithm': {'fusion_mode': True}},
    {'config': {'algorithm': {'fusion_mode': True}}},
])
def test_should_run_fusion_when_fusion_mode_and_default_caller(config):
    result = tophat._should_run_fusion(config)
    assert result is True


@pytest.mark.parametrize('config', [
    {'algorithm': {'fusion_mode': True, 'fusion_caller': 'TEST'}},
    {'config': {'algorithm': {'fusion_mode': True, 'fusion_caller': 'TEST'}}},
])
def test_should_not_run_fusion_when_fusion_mode_and_caller_not_tophat(config):
    result = tophat._should_run_fusion(config)
    assert result is False


@pytest.mark.parametrize('config', [
    {'algorithm': {'fusion_mode': True, 'fusion_caller': 'tophat'}},
    {'config': {'algorithm': {'fusion_mode': True, 'fusion_caller': 'tophat'}}},
])
def test_should_run_fusion_when_fusion_mode_and_caller_tophat(config):
    result = tophat._should_run_fusion(config)
    assert result is True


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
