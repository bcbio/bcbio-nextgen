from bcbio.pipeline import datadict as dd

import pytest


FUSION_CALLER = 'TEST_FUSION_CALLER'


@pytest.fixture
def data():
    return {
        'config': {
            'algorithm': {
                'fusion_caller': FUSION_CALLER,
            },
        },
    }


def test_get_fusion_caller(data):
    result = dd.get_fusion_caller(data)
    assert result == FUSION_CALLER
