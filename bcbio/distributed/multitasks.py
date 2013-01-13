"""Multiprocessing ready entry points for sample analysis.
"""
from bcbio import utils
from bcbio.pipeline import sample, lane, shared, variation
from bcbio.variation import realign, genotype, ensemble, recalibrate, multi

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
def prep_recal(*args):
    return recalibrate.prep_recal(*args)

@utils.map_wrap
def write_recal_bam(*args):
    return recalibrate.write_recal_bam(*args)

@utils.map_wrap
def realign_sample(*args):
    return realign.realign_sample(*args)

@utils.map_wrap
def split_variants_by_sample(*args):
    return multi.split_variants_by_sample(*args)

@utils.map_wrap
def postprocess_variants(*args):
    return sample.postprocess_variants(*args)

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
def variantcall_sample(*args):
    return genotype.variantcall_sample(*args)

@utils.map_wrap
def combine_variant_files(*args):
    return genotype.combine_variant_files(*args)

@utils.map_wrap
def detect_sv(*args):
    return variation.detect_sv(*args)

@utils.map_wrap
def combine_calls(*args):
    return ensemble.combine_calls(*args)
