import mock
import pytest

from bcbio.pipeline import config_utils


@pytest.mark.parametrize('config', [
    {'algorithm': {'fusion_mode': False}},
    {'config': {'algorithm': {'fusion_mode': False}}},
])
def test_should_not_run_fusion_when_fusion_mode_false(config):
    result = config_utils.should_run_fusion(mock.Mock(), config)
    assert result is False


@pytest.mark.parametrize('config', [
    {'algorithm': {'fusion_mode': True}},
    {'config': {'algorithm': {'fusion_mode': True}}},
])
def test_should_run_fusion_when_fusion_mode_and_default_caller(config):
    result = config_utils.should_run_fusion(mock.Mock(), config)
    assert result is True


@pytest.mark.parametrize('config', [
    {'algorithm': {'fusion_mode': True, 'fusion_caller': 'TEST'}},
    {'config': {'algorithm': {'fusion_mode': True, 'fusion_caller': 'TEST'}}},
])
def test_should_not_run_fusion_when_fusion_mode_and_wrong_caller(config):
    result = config_utils.should_run_fusion('OTHER', config)
    assert result is False


@pytest.mark.parametrize('config', [
    {'algorithm': {'fusion_mode': True, 'fusion_caller': 'TEST'}},
    {'config': {'algorithm': {'fusion_mode': True, 'fusion_caller': 'TEST'}}},
])
def test_should_run_fusion_when_fusion_mode_and_right_caller(config):
    result = config_utils.should_run_fusion('TEST', config)
    assert result is True
