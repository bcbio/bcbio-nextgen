import os
import functools

import pytest
import pandas as pd
from bcbio.rnaseq import ericscript
from tests.conftest import make_workdir


class ConfigCreator(object):

    INPUT_DATA_DIR = 'fusion/input'
    INPUT_FILENAMES = {
        'work_bam':  '1_2_Test1.trimmed.fq.gz',
        'fq_files': [
            '1_1_Test1.trimmed.fq.gz',
            '1_2_Test1.trimmed.fq.gz',
        ],
    }
    DISAMBIGUATE = {
        'config': {
            'algorithm': {
                'disambiguate': ['mm9'],
            },
        },
    }
    BASE_CONFIG = {
        'rgnames': {'lane': 'TEST_LANE'},
    }

    def get_config_with_disambiguate(self, data_dir=None, work_dir=None):
        config = self.get_config_without_disambiguate(
            data_dir=data_dir, work_dir=work_dir)
        config.update(self.DISAMBIGUATE)
        return config

    def get_config_without_disambiguate(self, data_dir=None, work_dir=None):
        join = self._get_join_datadir_fn(data_dir)
        filepaths_config = self._join_filepaths(join, work_dir)
        return self._merge_dicts(self.BASE_CONFIG, filepaths_config)

    def _get_join_datadir_fn(self, data_dir):
        input_data_dir = os.path.join(data_dir, os.pardir, self.INPUT_DATA_DIR)
        return functools.partial(os.path.join, input_data_dir)

    def _join_filepaths(self, join_fn, workdir):
        return {
            'work_bam': join_fn(self.INPUT_FILENAMES['work_bam']),
            'dirs': {'work': workdir},
            'files': map(join_fn, self.INPUT_FILENAMES['fq_files'])
        }

    def _merge_dicts(self, *args):
        result = {}
        for d in args:
            result.update(d)
        return result


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
    return pd.read_csv(fname, sep='\t').drop('EricScore',  axis=1)


def test_detect_fusions_with_ericscipt_without_disambiguate(
        install_test_files, data_dir):
    """Run an RNA-seq analysis and test fusion genes
    """
    with make_workdir() as work_dir:
        sample_config = ConfigCreator().get_config_without_disambiguate(
            data_dir=data_dir, work_dir=work_dir)
        ericscript.run(sample_config)
        assert_run_successfully(work_dir=work_dir, data_dir=data_dir)


def test_detect_fusions_with_ericscipt_with_disambiguate(
        install_test_files, data_dir):
    """Run an RNA-seq analysis and test fusion genes
    """
    with make_workdir() as work_dir:
        sample_config = ConfigCreator().get_config_with_disambiguate(
            data_dir=data_dir, work_dir=work_dir)
        ericscript.run(sample_config)
        assert_run_successfully(work_dir=work_dir, data_dir=data_dir)
