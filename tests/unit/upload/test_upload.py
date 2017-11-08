import os
import mock
import pytest

from bcbio import upload
from bcbio.rnaseq.ericscript import EricScriptConfig


@pytest.yield_fixture
def exists(mocker):
    yield mocker.patch.object(os.path, 'exists')


@pytest.yield_fixture
def es_out_dir(mocker):
    yield mocker.patch.object(EricScriptConfig, 'sample_out_dir')


def test_add_ericscript_files_does_nothing_if_oudir_doesnt_exist(
        exists, es_out_dir):
    exists.return_value = False
    result = upload._maybe_add_ericscript_files(mock.Mock(), mock.Mock(), [])
    assert result == []


def test_add_ericscript_files_appends_ES_data_to_out_if_outdir_exists(
        exists, es_out_dir):
    exists.return_value = True
    result = upload._maybe_add_ericscript_files(mock.Mock(), mock.Mock(), [])
    expected = [{
        'path': es_out_dir,
        'type': 'directory',
        'ext': 'ericscript',
    }]
    assert result == expected
