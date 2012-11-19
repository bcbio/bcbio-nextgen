"""Ipython parallel ready entry points for parallel execution
"""
from IPython.parallel import require

from bcbio.pipeline import sample, lane, shared, variation
from bcbio.variation import realign, genotype

@require(lane)
def process_lane(*args, **kwargs):
    return apply(lane.process_lane, *args, **kwargs)
