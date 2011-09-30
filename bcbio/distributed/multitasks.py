"""Multiprocessing ready entry points for sample analysis.
"""
from bcbio import utils
from bcbio.pipeline import sample, lane

@utils.map_wrap
def process_lane(*args):
    return lane.process_lane(*args)

@utils.map_wrap
def process_alignment(*args):
    return lane.process_alignment(*args)

@utils.map_wrap
def process_sample(*args):
    return sample.process_sample(*args)
