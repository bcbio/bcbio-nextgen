from copy import deepcopy
import functools
import os
import pytest

import pandas as pd

from bcbio.rnaseq import ericscript
from tests.conftest import make_workdir
from bcbio.pipeline import config_utils, run_info
from bcbio.log import setup_script_logging


def create_sample_config(data_dir, work_dir, disambiguate=False):
    system_config, system_file = config_utils.load_system_config(work_dir=work_dir)
    system_config["dirs"] = run_info.setup_directories(work_dir, work_dir, system_config, system_file)
    c = ConfigCreator(data_dir, work_dir, system_config)
    if disambiguate:
        return c.config_with_disambiguate()
    else:
        return c.config_without_disambiguate()


@pytest.fixture
def setup_logging():
    setup_script_logging()


@pytest.mark.ericscript
@pytest.mark.install_required
def test_detect_fusions_with_ericscipt_without_disambiguate(
        install_test_files, data_dir, setup_logging):
    """Run gene fusion analysis on trimmed pair-end reads with EricScript.
       Requires installation of EricScript and its reference data.
    """
    with make_workdir() as work_dir:
        sample_config = create_sample_config(
            data_dir, work_dir, disambiguate=False)
        ericscript.run(sample_config)
        assert_run_successfully(work_dir=work_dir, data_dir=data_dir)


@pytest.mark.ericscript
@pytest.mark.install_required
def test_detect_fusions_with_ericscipt_with_disambiguate(
        install_test_files, data_dir, setup_logging):
    """Run gene fusion analysis on disambiguated reads with EricScript.
       Requires installation of EricScript and its reference data.
    """
    with make_workdir() as work_dir:
        sample_config = create_sample_config(
            data_dir, work_dir, disambiguate=True)
        ericscript.run(sample_config)
        assert_run_successfully(work_dir=work_dir, data_dir=data_dir)


class ConfigCreator(object):

    _INPUT_DATA_DIR = 'fusion/input'
    _INPUT_FILENAMES = {
        'work_bam':  'Test1.nsorted.human.sorted.bam',
        'fq_files': [
            '1_1_Test1.trimmed.fq.gz',
            '1_2_Test1.trimmed.fq.gz',
        ],
    }

    def __init__(self, data_dir, work_dir, system_config):
        self._data_dir = data_dir
        self._work_dir = work_dir
        self._system_config = deepcopy(system_config)

    def config_without_disambiguate(self):
        config = self._get_base_config()
        config["genome_build"] = "hg19"
        config.update(self._get_filepaths())
        config = run_info.add_reference_resources(config)
        return config

    def config_with_disambiguate(self):
        config = self.config_without_disambiguate()
        if 'algorithm' not in config['config']:
            config['config']['algorithm'] = {}
        config['config']['algorithm'].update({'disambiguate': ['mm9']})
        return config

    def _get_base_config(self):
        conf = deepcopy(self._system_config)
        return {
            'analysis': 'rna-seq',
            'rgnames': {'lane': 'TEST_LANE'},
            'config': conf,
            'dirs': conf["dirs"],
        }

    def _get_filepaths(self):
        data_dir = os.path.join(
            self._data_dir, os.pardir, self._INPUT_DATA_DIR)
        join_fn = functools.partial(os.path.join, data_dir)
        return {
            'work_bam': join_fn(self._INPUT_FILENAMES['work_bam']),
            'files': map(join_fn, self._INPUT_FILENAMES['fq_files'])
        }


def assert_run_successfully(data_dir=None, work_dir=None):
    ERICSCRIPT_DIR = 'ericscript/TEST_LANE'
    OUT_DIR = 'out'
    ALN_DIR = 'aln'
    assert os.path.exists(os.path.join(work_dir, ERICSCRIPT_DIR))
    assert os.path.exists(os.path.join(ERICSCRIPT_DIR, OUT_DIR))
    assert os.path.exists(os.path.join(ERICSCRIPT_DIR, ALN_DIR))

    result = os.path.join(
        work_dir, ERICSCRIPT_DIR, 'TEST_LANE.results.total.tsv')
    expected_result = os.path.join(
        data_dir, os.pardir, 'fusion/results/TEST_LANE.results.total.tsv')
    expected = _load_result_file(expected_result)
    result = _load_result_file(result)
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
    df = pd.read_csv(fname, sep='\t')[columns_to_keep].sort_values(by=sort_by)
    df.index = range(len(df))
    return df
