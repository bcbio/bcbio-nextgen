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


def test_detect_fusion_callers(ericscript_run, sample):
    samples = [[sample], [sample]]
    rnaseq.detect_fusions(samples)
    assert ericscript_run.call_count == len(samples)


