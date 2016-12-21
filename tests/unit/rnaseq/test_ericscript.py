import pytest
import mock

from bcbio.rnaseq import ericscript
from bcbio.rnaseq.ericscript import EricScriptConfig
from tests.unit.conftest import DummyFileTransaction


@pytest.yield_fixture
def utils(mocker):
    yield mocker.patch('bcbio.rnaseq.ericscript.utils')


class TestEricScriptConfig(object):

    @pytest.yield_fixture
    def es_config(self, utils):
        sample = {
            'rgnames': {'lane': 'TEST_LANE'},
            'dirs': {
                'work': 'TEST_WORK_DIR'
            },
            'config': {
                'resources': {'ericscript': {
                    'env': '/path/to/envs/ericscript',
                    'db': '/path/to/ericscript_db',
                    },
                }
            }
        }
        yield EricScriptConfig(sample)

    def test_info_message(self, es_config):
        assert es_config.info_message == "Detect gene fusions with EricScript"

    def test_gets_env_via_utils(self, es_config):
        expected_env = ericscript.utils.get_ericscript_env.return_value
        assert es_config.env == expected_env

    def test_calls_utils_with_env_prefix_to_get_env(self, es_config):
        es_config.env
        ericscript.utils.get_ericscript_env.assert_called_once_with(
            '/path/to/envs/ericscript')

    def test_output_dir(self, es_config):
        assert es_config.output_dir == 'TEST_WORK_DIR/ericscript'

    def test_sample_output_dir(self, es_config):
        assert es_config.sample_out_dir == 'TEST_WORK_DIR/ericscript/TEST_LANE'

    def test_reference_index(self, es_config):
        expected = '/path/to/ericscript_db/data/homo_sapiens/allseq.fa.bwt'
        assert es_config.reference_index == expected

    def test_reference_fasta(self, es_config):
        expected = '/path/to/ericscript_db/data/homo_sapiens/allseq.fa'
        assert es_config.reference_fasta == expected

    def test_get_run_command(self, es_config):
        tx_dir = 'TX_DIR'
        input_files = ('file1.fq', 'file2.fq')
        cmd = es_config.get_run_command(tx_dir, input_files)
        expected = [
            'ericscript.pl',
            '-db',
            '/path/to/ericscript_db',
            '-name',
            'TEST_LANE',
            '-o',
            'TX_DIR',
            'file1.fq',
            'file2.fq',
        ]
        assert cmd == expected


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
    def es_config(self, mocker):
        mock_ES = mocker.patch(
            'bcbio.rnaseq.ericscript.EricScriptConfig',
            autospec=True
        )
        yield mock_ES(mock.Mock())

    @pytest.yield_fixture
    def do_run(self, mocker):
        yield mocker.patch(
            'bcbio.rnaseq.ericscript.do.run',
            autospec=True
        )

    @pytest.yield_fixture
    def prepare_data(self, mocker):
        yield mocker.patch('bcbio.rnaseq.ericscript.prepare_input_data')

    def test_returns_sample_config(
            self, prepare_data, mock_ft, do_run, es_config, utils):
        config = mock.MagicMock()
        result = ericscript.run(config)
        assert result == config

    def test_gets_ericscript_command(
            self, prepare_data, mock_ft, do_run, es_config, utils):

        ericscript.run(mock.Mock())
        es_config.get_run_command.assert_called_once_with(
            mock.ANY,
            prepare_data.return_value
        )

    def test_runs_ericscript_command(
            self, prepare_data, mock_ft, do_run, es_config, utils):
        ericscript.run(mock.Mock())
        do_run.assert_called_once_with(
            es_config.get_run_command.return_value,
            es_config.info_message,
            env=es_config.env,
        )

    def test_calls_file_transaction(
            self, prepare_data, mock_ft, do_run, es_config, utils):
        config = mock.MagicMock()
        ericscript.run(config)
        mock_ft.assert_called_once_with(config, es_config.sample_out_dir)

    def test_creates_base_ericscript_output_dir(
            self, prepare_data, mock_ft, do_run, es_config, utils):
        config = mock.MagicMock()
        ericscript.run(config)
        utils.safe_makedir.assert_called_once_with(es_config.output_dir)


class TestBuildBWAIndex(object):

    @pytest.fixture
    def exists(self, mocker):
        yield mocker.patch('bcbio.rnaseq.ericscript.os.path.exists')

    @pytest.fixture
    def bwa_build(self, mocker):
        yield mocker.patch('bcbio.rnaseq.ericscript.bwa.build_bwa_index')

    def test_doesnt_build_index_if_already_exists(self, exists, bwa_build):
        exists.return_value = True
        ericscript.build_bwa_index_if_absent(mock.Mock(), mock.Mock())
        assert bwa_build.call_count == 0

    def test_builds_index_if_it_doesnt_exist(self, exists, bwa_build):
        exists.return_value = False
        es_config, config = mock.Mock(spec=EricScriptConfig), mock.Mock()
        ericscript.build_bwa_index_if_absent(es_config, config)
        bwa_build.assert_called_once_with(es_config.reference_fasta, config)

    def test_checks_if_index_exists(self, exists, bwa_build):
        es_config = mock.Mock(spec=EricScriptConfig)
        ericscript.build_bwa_index_if_absent(es_config, mock.Mock())
        exists.assert_called_once_with(es_config.reference_index)
