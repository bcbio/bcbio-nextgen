import pytest

from bcbio.rnaseq import ericscript
from tests.unit.conftest import DummyFileTransaction


@pytest.yield_fixture
def utils(mocker):
    yield mocker.patch('bcbio.rnaseq.ericscript.utils')

class TestGetInputData(object):
    def test_get_disambiguated_bam(self, mocker):
        sample_config = {
            'config': {
                'algorithm': {
                    'disambiguate': ['mm9'],
                }
            },
            'work_bam':
                '/path/to/disambiguate_star/Test1.nsorted.human.sorted.bam',
            'dirs': {'work': '/path/to/workdir'},
        }

        convert = mocker.patch('bcbio.rnaseq.ericscript.convert_bam_to_fastq')
        result = ericscript.prepare_input_data(sample_config)
        convert.assert_called_once_with(
            sample_config['work_bam'],
            sample_config['dirs']['work'],
            None, None, sample_config
        )
        assert result == convert.return_value

    def test_get_fastq_input_files_if_no_disambiguation(self):
        fq_files = (
            '/path/to/1_1_trimmed.fq.gz',
            '/path/to/1_2_trimmed.fq.gz'
        )
        sample_config = {'files': list(fq_files)}
        result = ericscript.prepare_input_data(sample_config)
        assert result == fq_files


class TestRun(object):

    @pytest.yield_fixture
    def mock_ft(self, mocker):
        yield mocker.patch(
            'bcbio.rnaseq.ericscript.file_transaction',
            side_effect=DummyFileTransaction
        )

    @pytest.yield_fixture
    def do_run(self, mocker):
        yield mocker.patch(
            'bcbio.rnaseq.ericscript.do.run',
            autospec=True
        )

    @pytest.yield_fixture
    def prepare_data(self, mocker):
        yield mocker.patch('bcbio.rnaseq.ericscript.prepare_input_data')
