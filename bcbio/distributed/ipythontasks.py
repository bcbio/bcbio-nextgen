"""Ipython parallel ready entry points for parallel execution
"""
import contextlib
import json
import zlib

from IPython.parallel import require

from bcbio import chipseq, structural
from bcbio.bam import callable
from bcbio.ngsalign import alignprep
from bcbio.pipeline import (archive, config_utils, disambiguate, sample, lane, qcsummary, shared,
                            variation, rnaseq)
from bcbio.provenance import system
from bcbio.variation import (bamprep, coverage, genotype, ensemble, multi, population,
                             recalibrate, validate, vcfutils)
from bcbio.log import logger, setup_local_logging

@contextlib.contextmanager
def _setup_logging(args):
    config = None
    if len(args) == 1 and isinstance(args[0], (list, tuple)):
        args = args[0]
    for arg in args:
        if config_utils.is_nested_config_arg(arg):
            config = arg["config"]
            break
        elif config_utils.is_std_config_arg(arg):
            config = arg
            break
        elif isinstance(arg, (list, tuple)) and config_utils.is_nested_config_arg(arg[0]):
            config = arg[0]["config"]
            break
    if config is None:
        raise NotImplementedError("No config found in arguments: %s" % args[0])
    handler = setup_local_logging(config, config.get("parallel", {}))
    try:
        yield None
    except:
        logger.exception("Unexpected error")
        raise
    finally:
        if hasattr(handler, "close"):
            handler.close()

def _unzip_args(args):
    """Unzip arguments if passed as compressed JSON string.
    """
    if len(args) == 1 and isinstance(args[0], basestring):
        return [json.loads(zlib.decompress(args[0]))]
    else:
        return args

@require(lane)
def process_lane(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(lane.process_lane, *args)

@require(lane)
def trim_lane(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(lane.trim_lane, *args)

@require(lane)
def process_alignment(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(lane.process_alignment, *args)

@require(alignprep)
def prep_align_inputs(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(alignprep.create_inputs, *args)

@require(lane)
def postprocess_alignment(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(lane.postprocess_alignment, *args)

@require(sample)
def merge_sample(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(sample.merge_sample, *args)

@require(sample)
def delayed_bam_merge(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(sample.delayed_bam_merge, *args)

@require(sample)
def recalibrate_sample(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(sample.recalibrate_sample, *args)

@require(recalibrate)
def prep_recal(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(recalibrate.prep_recal, *args)

@require(multi)
def split_variants_by_sample(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(multi.split_variants_by_sample, *args)

@require(bamprep)
def piped_bamprep(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(bamprep.piped_bamprep, *args)

@require(variation)
def postprocess_variants(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(variation.postprocess_variants, *args)

@require(qcsummary)
def pipeline_summary(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(qcsummary.pipeline_summary, *args)

@require(rnaseq)
def generate_transcript_counts(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(rnaseq.generate_transcript_counts, *args)

@require(rnaseq)
def run_cufflinks(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(rnaseq.run_cufflinks, *args)

@require(shared)
def combine_bam(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(shared.combine_bam, *args)

@require(callable)
def combine_sample_regions(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(callable.combine_sample_regions, *args)

@require(genotype)
def variantcall_sample(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(genotype.variantcall_sample, *args)

@require(vcfutils)
def combine_variant_files(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(vcfutils.combine_variant_files, *args)

@require(vcfutils)
def concat_variant_files(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(vcfutils.concat_variant_files, *args)

@require(vcfutils)
def merge_variant_files(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(vcfutils.merge_variant_files, *args)

@require(population)
def prep_gemini_db(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(population.prep_gemini_db, *args)

@require(structural)
def detect_sv(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(structural.detect_sv, *args)

@require(ensemble)
def combine_calls(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(ensemble.combine_calls, *args)

@require(validate)
def compare_to_rm(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(validate.compare_to_rm, *args)

@require(coverage)
def coverage_summary(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(coverage.summary, *args)

@require(disambiguate)
def run_disambiguate(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(disambiguate.run, *args)

@require(system)
def machine_info(*args):
    args = _unzip_args(args)
    return system.machine_info()

@require(chipseq)
def clean_chipseq_alignment(*args):
    args = _unzip_args(args)
    return chipseq.machine_info()

@require(archive)
def archive_to_cram(*args):
    args = _unzip_args(args)
    with _setup_logging(args):
        return apply(archive.to_cram, *args)
