from bcbio.ngsalign import tophat


def test_set_fusion_mode_false():
    config = {
        'algorithm': {
            'fusion_mode': False
        }
    }
    opts = {}
    tophat._set_fusion_mode(opts, config)
    assert opts == {}


def test_set_fusion_mode_true__default_caller():
    config = {
        'algorithm': {
            'fusion_mode': True
        }
    }
    opts = {}
    tophat._set_fusion_mode(opts, config)
    assert opts == {'fusion-search': True}


def test_set_fusion_mode_true__wrong_caller():
    config = {
        'algorithm': {
            'fusion_mode': True,
            'fusion_caller': 'ericscript',
        }
    }
    opts = {}
    tophat._set_fusion_mode(opts, config)
    assert opts == {}


def test_set_fusion_mode_true__tophat_caller():
    config = {
        'algorithm': {
            'fusion_mode': True,
            'fusion_caller': 'tophat',
        }
    }
    opts = {}
    tophat._set_fusion_mode(opts, config)
    assert opts == {'fusion-search': True}
