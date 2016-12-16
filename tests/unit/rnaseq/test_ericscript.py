import pytest
import mock

from bcbio.rnaseq import ericscript
from bcbio.rnaseq.ericscript import EricScriptConfig
from tests.unit.conftest import DummyFileTransaction

class TestEricScriptConfig(object):

    @pytest.yield_fixture
    def utils(self, mocker):
        yield mocker.patch('bcbio.rnaseq.ericscript.utils')

    @pytest.yield_fixture
    def es_config(self, utils):
        sample = {
            'rgnames': {'lane': 'TEST_LANE'},
            'dirs': {
                'work': 'TEST_WORK_DIR'
            },
        }
        yield EricScriptConfig(sample)

    def test_info_message(self, es_config):
        assert es_config.info_message == "Detect gene fusions with EricScript"

    def test_env(self, es_config):
        expected_env = ericscript.utils.get_ericscript_env.return_value
        assert es_config.env == expected_env

    def test_output_dir(self, es_config):
        assert es_config.output_dir == 'TEST_WORK_DIR/ericscript'

    def test_get_run_command(self, es_config):
        tx_dir = 'TX_DIR'
        input_files = ('file1.fq', 'file2.fq')
        cmd = es_config.get_run_command(tx_dir, input_files)
        expected = [
            'ericscript.pl',
            '-db',
            '/data/ericscript/ericscript_db',
            '-name',
            'TEST_LANE',
            '-o',
            'TX_DIR',
            'file1.fq',
            'file2.fq',
        ]
        assert cmd == expected


class TestGetInputData(object):
    def test_get_disambiguated_bam(self):
        sample_config = {
            'config': {
                'algorithm': {
                    'disambiguate': ['mm9'],
                }
            },
            'work_bam':
                '/path/to/disambiguate_star/Test1.nsorted.human.sorted.bam',
        }

        result = ericscript.get_input_data('TX_DIR', sample_config)
        expected_fq_files = ['TX_DIR/input1.fq', 'TX_DIR/input2.fq']
        assert result.files == expected_fq_files
        assert result.convert_cmd == [
            'bamToFastq',
            '-i',
            '/path/to/disambiguate_star/Test1.nsorted.human.sorted.bam',
            '-fq',
            expected_fq_files[0],
            '-fq2',
            expected_fq_files[1],
        ]

    def test_get_fastq_input_files_if_no_disambiguation(self):
        fq_files = (
            '/path/to/1_1_trimmed.fq.gz',
            '/path/to/1_2_trimmed.fq.gz'
        )
        sample_config = {'files': list(fq_files)}
        result = ericscript.get_input_data('TX_DIR', sample_config)
        assert result.files == fq_files
        assert result.convert_cmd is None

class TestRun(object):

    @pytest.yield_fixture
    def mock_ft(self, mocker):
        yield mocker.patch(
            'bcbio.rnaseq.ericscript.file_transaction',
            side_effect=DummyFileTransaction
        )

    @pytest.yield_fixture
    def es_config(self, mocker):
        mock_ES = mocker.patch(
            'bcbio.rnaseq.ericscript.EricScriptConfig', autospec=True)
        yield mock_ES(mock.Mock())

    @pytest.yield_fixture
    def do_run(self, mocker):
        yield mocker.patch(
            'bcbio.rnaseq.ericscript.do.run', autospec=True)

    def test_returns_sample_config(self, mock_ft, do_run, es_config):
        config = mock.MagicMock()
        result = ericscript.run(config)
        assert result == config

    def test_run_ericscript_without_input_file_conversion(
            self, mocker, do_run, es_config, mock_ft):
        get_data = mocker.patch('bcbio.rnaseq.ericscript.get_input_data')
        get_data.return_value = ericscript.InputData('FILES', None)

        ericscript.run(mock.Mock())
        do_run.assert_called_once_with(
            es_config.get_run_command.return_value,
            es_config.info_message,
            env=es_config.env,
        )

    def test_run_ericscript_with_input_file_conversion(
            self, mocker, do_run, es_config, mock_ft):
        get_data = mocker.patch('bcbio.rnaseq.ericscript.get_input_data')
        get_data.return_value = ericscript.InputData('FILES', 'CONVERT_CMD')

        ericscript.run(mock.Mock())
        # assert that both input file conversion and EricScript were run:
        convert_call = mock.call(
            'CONVERT_CMD',
            'Convert disambiguated bam reads to fastq',
            env=es_config.env,
        )
        ericscript_call = mock.call(
            es_config.get_run_command.return_value,
            es_config.info_message,
            env=es_config.env,
        )
        do_run.assert_has_calls([convert_call, ericscript_call])

    def test_calls_file_transaction(self, do_run, mock_ft, es_config):
        config = mock.MagicMock()
        ericscript.run(config)
        mock_ft.assert_called_once_with(config, es_config.output_dir)
