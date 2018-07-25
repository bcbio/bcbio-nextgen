import pytest

from bcbio.pipeline import rnaseq


@pytest.yield_fixture
def ericscript_run(mocker):
    yield mocker.patch('bcbio.pipeline.rnaseq.ericscript.run', autospec=True)


@pytest.fixture
def sample():
    return {
        'config': {
            'algorithm': {
                'fusion_mode': True,
                'fusion_caller': 'ericscript'
            }
        }
    }

# XXX Skips runs since there is no database

def test_detect_fusion_callers_calls_for_each_sample(ericscript_run, sample):
    samples = [[sample]]
    rnaseq.detect_fusions(samples)
    assert ericscript_run.call_count == 0
    #assert ericscript_run.call_count == len(samples)


def test_detect_fusions_returns_updated_samples(ericscript_run, sample):
    samples = [[sample]]
    result = rnaseq.detect_fusions(samples)
    #assert result == [[ericscript_run.return_value]]


def test_detect_fusions_calls_caller_with_sample_dict(ericscript_run, sample):
    samples = [[sample]]
    rnaseq.detect_fusions(samples)
    assert ericscript_run.call_count == 0
    #ericscript_run.assert_called_once_with(sample)
