"""Ipython parallel ready entry points for parallel execution
"""
import contextlib

from IPython.parallel import require

from bcbio.distributed import ipython
from bcbio.pipeline import sample, lane, qcsummary, shared, variation
from bcbio.provenance import system
from bcbio.variation import (bamprep, coverage, realign, genotype, ensemble, multi, population,
                             recalibrate, validate, vcfutils)
from bcbio.log import logger, setup_local_logging

@contextlib.contextmanager
def _setup_logging(args):
    config = None
    if len(args) == 1 and isinstance(args[0], (list, tuple)):
        args = args[0]
    for arg in args:
        if ipython.is_nested_config_arg(arg):
            config = arg["config"]
            break
        elif ipython.is_std_config_arg(arg):
            config = arg
            break
    if config is None:
        raise NotImplementedError("No config in %s:" % args[0])
    handler = setup_local_logging(config, config.get("parallel", {}))
    try:
        yield None
        if hasattr(handler, "close"):
            handler.close()
    except:
        logger.exception("Unexpected error")
        raise

@require(lane)
def process_lane(*args):
    with _setup_logging(args):
        return apply(lane.process_lane, *args)

@require(lane)
def trim_lane(*args):
    with _setup_logging(args):
        return apply(lane.trim_lane, *args)

@require(lane)
def process_alignment(*args):
    with _setup_logging(args):
        return apply(lane.process_alignment, *args)
process_alignment.metadata = {"resources": ["novoalign", "bwa", "bowtie", "tophat"]}

@require(lane)
def align_prep_full(*args):
    with _setup_logging(args):
        return apply(lane.align_prep_full, *args)
align_prep_full.metadata = {"resources": ["novoalign", "bwa", "gatk"]}

@require(sample)
def merge_sample(*args):
    with _setup_logging(args):
        return apply(sample.merge_sample, *args)

@require(sample)
def delayed_bam_merge(*args):
    with _setup_logging(args):
        return apply(sample.delayed_bam_merge, *args)
delayed_bam_merge.metadata = {"resources": ["samtools"]}

@require(sample)
def recalibrate_sample(*args):
    with _setup_logging(args):
        return apply(sample.recalibrate_sample, *args)

@require(recalibrate)
def prep_recal(*args):
    with _setup_logging(args):
        return apply(recalibrate.prep_recal, *args)
prep_recal.metadata = {"resources": ["gatk"]}

@require(recalibrate)
def write_recal_bam(*args):
    with _setup_logging(args):
        return apply(recalibrate.write_recal_bam, *args)

@require(realign)
def realign_sample(*args):
    with _setup_logging(args):
        return apply(realign.realign_sample, *args)

@require(multi)
def split_variants_by_sample(*args):
    with _setup_logging(args):
        return apply(multi.split_variants_by_sample, *args)

@require(bamprep)
def piped_bamprep(*args):
    with _setup_logging(args):
        return apply(bamprep.piped_bamprep, *args)

@require(variation)
def postprocess_variants(*args):
    with _setup_logging(args):
        return apply(variation.postprocess_variants, *args)
postprocess_variants.metadata = {"resources": ["gatk-vqsr", "gatk", "snpEff"]}

@require(qcsummary)
def pipeline_summary(*args):
    with _setup_logging(args):
        return apply(qcsummary.pipeline_summary, *args)

@require(sample)
def generate_transcript_counts(*args):
    with _setup_logging(args):
        return apply(sample.generate_transcript_counts, *args)

@require(sample)
def generate_bigwig(*args):
    with _setup_logging(args):
        return apply(sample.generate_bigwig, *args)

@require(shared)
def combine_bam(*args):
    with _setup_logging(args):
        return apply(shared.combine_bam, *args)

@require(genotype)
def variantcall_sample(*args):
    with _setup_logging(args):
        return apply(genotype.variantcall_sample, *args)
variantcall_sample.metadata = {"resources": ["gatk", "freebayes", "gatk-haplotype"]}

@require(vcfutils)
def combine_variant_files(*args):
    with _setup_logging(args):
        return apply(vcfutils.combine_variant_files, *args)

@require(vcfutils)
def concat_variant_files(*args):
    with _setup_logging(args):
        return apply(vcfutils.concat_variant_files, *args)

@require(population)
def prep_gemini_db(*args):
    with _setup_logging(args):
        return apply(population.prep_gemini_db, *args)
prep_gemini_db.metadata = {"resources": ["gemini"]}

@require(variation)
def detect_sv(*args):
    with _setup_logging(args):
        return apply(variation.detect_sv, *args)

@require(ensemble)
def combine_calls(*args):
    with _setup_logging(args):
        return apply(ensemble.combine_calls, *args)

@require(validate)
def compare_to_rm(*args):
    with _setup_logging(args):
        return apply(validate.compare_to_rm, *args)

@require(coverage)
def coverage_summary(*args):
    with _setup_logging(args):
        return apply(coverage.summary, *args)
coverage_summary.metadata = {"resources": ["bcbio_coverage"]}

@require(system)
def machine_info(*args):
    return system.machine_info()
