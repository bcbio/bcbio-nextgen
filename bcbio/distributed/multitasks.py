"""Multiprocessing ready entry points for sample analysis.
"""
from bcbio import utils
from bcbio.pipeline import sample, lane, shared
from bcbio.variation import realign, genotype

@utils.map_wrap
def process_lane(*args):
    return lane.process_lane(*args)

@utils.map_wrap
def process_alignment(*args):
    return lane.process_alignment(*args)

@utils.map_wrap
def merge_sample(*args):
    return sample.merge_sample(*args)

@utils.map_wrap
def recalibrate_sample(*args):
    return sample.recalibrate_sample(*args)

@utils.map_wrap
def realign_sample(*args):
    return realign.realign_sample(*args)

@utils.map_wrap
def process_sample(*args):
    return sample.process_sample(*args)

@utils.map_wrap
def generate_bigwig(*args):
    return sample.generate_bigwig(*args)

@utils.map_wrap
def combine_bam(*args):
    return shared.combine_bam(*args)

@utils.map_wrap
def unified_genotyper_sample(*args):
    return genotype.unified_genotyper_sample(*args)

@utils.map_wrap
def combine_variant_files(*args):
    return genotype.combine_variant_files(*args)
