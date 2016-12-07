from tests.unit.data import DATA as _DATA
from tests.unit.data import NAMES as _NAMES
from copy import deepcopy

import pytest

from bcbio.ngsalign.star import _get_star_dirnames
from bcbio.pipeline import datadict as dd


@pytest.fixture
def data():
    return deepcopy(_DATA)


@pytest.fixture
def names():
    return deepcopy(_NAMES)


def test_get_star_dirnames(data, names):
    align_dir = '/path/to/align/dir'
    lane = dd.get_lane(data)
    result = _get_star_dirnames(align_dir, data, names)
    assert result.out_dir == '/path/to/align/dir/%s_star' % lane
    assert result.out_prefix == '/path/to/align/dir/%s' % lane
    assert result.out_file == '/path/to/align/dir/%sAligned.out.sam' % lane
    assert result.final_out == '/path/to/align/dir/%s_star/%s.bam' % (
        lane, names['sample'])
