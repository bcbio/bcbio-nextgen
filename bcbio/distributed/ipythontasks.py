"""Ipython parallel ready entry points for parallel execution
"""
from IPython.parallel import require

from bcbio.pipeline import sample, lane, shared, variation
from bcbio.variation import realign, genotype

@require(lane)
def process_lane(*args, **kwargs):
    return apply(lane.process_lane, *args, **kwargs)

@require(lane)
def process_alignment(*args):
    return apply(lane.process_alignment, *args)

@require(sample)
def merge_sample(*args):
    return apply(sample.merge_sample, *args)

@require(sample)
def recalibrate_sample(*args):
    return apply(sample.recalibrate_sample, *args)

@require(realign)
def realign_sample(*args):
    return apply(realign.realign_sample, *args)

@require(sample)
def postprocess_variants(*args):
    return apply(sample.postprocess_variants, *args)

@require(sample)
def process_sample(*args):
    return apply(sample.process_sample, *args)

@require(sample)
def generate_bigwig(*args):
    return apply(sample.generate_bigwig, *args)

@require(shared)
def combine_bam(*args):
    return apply(shared.combine_bam, *args)

@require(genotype)
def variantcall_sample(*args):
    return apply(genotype.variantcall_sample, *args)

@require(genotype)
def combine_variant_files(*args):
    return apply(genotype.combine_variant_files, *args)

@require(variation)
def detect_sv(*args):
    return apply(variation.detect_sv, *args)
