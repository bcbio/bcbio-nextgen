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


def test_get_ericscript_env():
    sample = {
        'config': {
            'resources': {'ericscript': {
                'env': '/path/to/envs/ericscript'
                },
            }
        }
    }
    result = dd.get_ericscript_env(sample)
    assert result == '/path/to/envs/ericscript'


def test_get_ericscript_db():
    sample = {
        'config': {
            'resources': {'ericscript': {
                'db': '/path/to/ericscript_db'
                },
            }
        }
    }
    result = dd.get_ericscript_db(sample)
    assert result == '/path/to/ericscript_db'
