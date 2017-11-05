from bcbio.pipeline import datadict as dd


def test_get_fusion_caller():
    data = {
        'config': {
            'algorithm': {
                'fusion_caller': 'FUSION_CALLER',
            },
        },
    }

    result = dd.get_fusion_caller(data)
    assert result == 'FUSION_CALLER'
