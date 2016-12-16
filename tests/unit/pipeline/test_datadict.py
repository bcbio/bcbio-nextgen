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


def test_get_ericscript_outdir():
    sample_config = {
        'dirs': {
            'work': 'TEST_WORK_DIR'
        },
    }
    result = dd.get_ericscript_outdir(sample_config)
    expected = 'TEST_WORK_DIR/ericscript'
    assert result == expected
