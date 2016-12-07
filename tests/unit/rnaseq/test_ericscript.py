from copy import deepcopy

import pytest

from bcbio.rnaseq import ericscript
from tests.unit.data import CONFIG as _CONFIG


@pytest.fixture
def config():
    return deepcopy(_CONFIG[0][0])


@pytest.fixture
def mock_io(mocker):
    mocker.patch('bcbio.rnaseq.ericscript.file_transaction')


def test_get_fastq_files(config):
    result = ericscript._get_fastq_files(config)
    expected = (
        '/bcbio-nextgen/tests/test_automated_output/trimmed/1_1_Test1.trimmed.fq.gz',  # noqa
        '/bcbio-nextgen/tests/test_automated_output/trimmed/1_2_Test1.trimmed.fq.gz'   # noqa
    )

    assert result == expected
