import functools
import os

import pandas as pd
from bcbio.rnaseq import ericscript
from tests.conftest import make_workdir


class ConfigCreator(object):

    _INPUT_DATA_DIR = 'fusion/input'
    _INPUT_FILENAMES = {
        'work_bam':  'Test1.nsorted.human.sorted.bam',
        'fq_files': [
            '1_1_Test1.trimmed.fq.gz',
            '1_2_Test1.trimmed.fq.gz',
        ],
    }

    def get_config_with_disambiguate(self, data_dir=None, work_dir=None):
        config = self.get_config_without_disambiguate(
            data_dir=data_dir, work_dir=work_dir)
        config.update(self._get_disambiguate())
        return config

    def get_config_without_disambiguate(self, data_dir=None, work_dir=None):
        config = self._get_base_config()
        config.update(self._get_filepaths(data_dir, work_dir))
        return config

    def _get_base_config(self):
        return {
            'rgnames': {'lane': 'TEST_LANE'},
        }

    def _get_disambiguate(self):
        return {
            'config': {
                'algorithm': {'disambiguate': ['mm9']},
            },
        }

    def _get_filepaths(self, data_dir, workdir):
        data_dir = os.path.join(data_dir, os.pardir, self._INPUT_DATA_DIR)
        join_fn = functools.partial(os.path.join, data_dir)
        return {
            'work_bam': join_fn(self._INPUT_FILENAMES['work_bam']),
            'dirs': {'work': workdir},
            'files': map(join_fn, self._INPUT_FILENAMES['fq_files'])
        }


def assert_run_successfully(data_dir=None, work_dir=None):
    ERICSCRIPT_DIR = os.path.join(work_dir, 'ericscript')
    OUT_DIR = os.path.join(ERICSCRIPT_DIR, 'out')
    ALN_DIR = os.path.join(ERICSCRIPT_DIR, 'aln')
    EXPECTED_RESULT = os.path.join(
        data_dir, os.pardir, 'fusion/results/TEST_LANE.results.total.tsv')
    RESULT = os.path.join(work_dir, 'ericscript/TEST_LANE.results.total.tsv')

    assert os.path.exists(ERICSCRIPT_DIR)
    assert os.path.exists(OUT_DIR)
    assert os.path.exists(ALN_DIR)
    assert os.path.exists(RESULT)

    expected = _load_result_file(EXPECTED_RESULT)
    result = _load_result_file(RESULT)
    assert result.equals(expected)


def _load_result_file(fname):
    # Load results from the fname and check that the same fusions were
    # detected.
    # To compare dataframes, row order must be the same, so we sort values
    # and reset index.
    columns_to_keep = [
        'chr1',
        'chr2',
        'strand1',
        'strand2',
        'EnsemblGene1',
        'EnsemblGene2',
        'fusiontype'
    ]
    sort_by = ['EnsemblGene1', 'EnsemblGene2']
    df = pd.read_csv(fname, sep='\t')[columns_to_keep].sort(sort_by)
    df.index = range(len(df))
    return df


def test_detect_fusions_with_ericscipt_without_disambiguate(
        install_test_files, data_dir):
    """Run gene fusion analysis on trimmed pair-end reads with EricScript.
    """
    with make_workdir() as work_dir:
        sample_config = ConfigCreator().get_config_without_disambiguate(
            data_dir=data_dir, work_dir=work_dir)
        ericscript.run(sample_config)
        assert_run_successfully(work_dir=work_dir, data_dir=data_dir)


def test_detect_fusions_with_ericscipt_with_disambiguate(
        install_test_files, data_dir):
    """Run gene fusion analysis on disambiguated reads with EricScript.
    """
    with make_workdir() as work_dir:
        sample_config = ConfigCreator().get_config_with_disambiguate(
            data_dir=data_dir, work_dir=work_dir)
        ericscript.run(sample_config)
        assert_run_successfully(work_dir=work_dir, data_dir=data_dir)
