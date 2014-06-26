"""Ipython parallel ready entry points for parallel execution
"""
import contextlib

from IPython.parallel import require

from bcbio import chipseq, structural
from bcbio.bam import callable
from bcbio.rnaseq import sailfish
from bcbio.distributed import ipython
from bcbio.ngsalign import alignprep
from bcbio.pipeline import (archive, config_utils, disambiguate, sample, lane, qcsummary, shared,
                            variation, rnaseq)
from bcbio.provenance import system
from bcbio.variation import (bamprep, coverage, genotype, ensemble, joint, multi, population,
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
        yield config
    except:
        logger.exception("Unexpected error")
        raise
    finally:
        if hasattr(handler, "close"):
            handler.close()

@require(lane)
def process_lane(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(lane.process_lane, *args), config)

@require(lane)
def trim_lane(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(lane.trim_lane, *args), config)

@require(sailfish)
def run_sailfish(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args):
        return apply(sailfish.run_sailfish, *args)

@require(lane)
def process_alignment(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(lane.process_alignment, *args), config)

@require(alignprep)
def prep_align_inputs(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(alignprep.create_inputs, *args), config)

@require(lane)
def postprocess_alignment(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(lane.postprocess_alignment, *args), config)

@require(sample)
def merge_sample(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(sample.merge_sample, *args), config)

@require(sample)
def delayed_bam_merge(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(sample.delayed_bam_merge, *args), config)

@require(sample)
def recalibrate_sample(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(sample.recalibrate_sample, *args), config)

@require(recalibrate)
def prep_recal(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(recalibrate.prep_recal, *args), config)

@require(multi)
def split_variants_by_sample(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(multi.split_variants_by_sample, *args), config)

@require(bamprep)
def piped_bamprep(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(bamprep.piped_bamprep, *args), config)

@require(variation)
def postprocess_variants(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(variation.postprocess_variants, *args), config)

@require(qcsummary)
def pipeline_summary(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(qcsummary.pipeline_summary, *args), config)

@require(rnaseq)
def generate_transcript_counts(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(rnaseq.generate_transcript_counts, *args), config)

@require(rnaseq)
def run_cufflinks(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(rnaseq.run_cufflinks, *args), config)

@require(shared)
def combine_bam(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(shared.combine_bam, *args), config)

@require(callable)
def combine_sample_regions(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(callable.combine_sample_regions, *args), config)

@require(genotype)
def variantcall_sample(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(genotype.variantcall_sample, *args), config)

@require(vcfutils)
def combine_variant_files(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(vcfutils.combine_variant_files, *args), config)

@require(vcfutils)
def concat_variant_files(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(vcfutils.concat_variant_files, *args), config)

@require(vcfutils)
def merge_variant_files(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(vcfutils.merge_variant_files, *args), config)

@require(population)
def prep_gemini_db(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(population.prep_gemini_db, *args), config)

@require(structural)
def detect_sv(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(structural.detect_sv, *args), config)

@require(ensemble)
def combine_calls(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(ensemble.combine_calls, *args), config)

@require(validate)
def compare_to_rm(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(validate.compare_to_rm, *args), config)

@require(coverage)
def coverage_summary(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(coverage.summary, *args), config)

@require(disambiguate)
def run_disambiguate(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(disambiguate.run, *args), config)

@require(system)
def machine_info(*args):
    args = ipython.unzip_args(args)
    return system.machine_info()

@require(chipseq)
def clean_chipseq_alignment(*args):
    args = ipython.unzip_args(args)
    return chipseq.machine_info()

@require(archive)
def archive_to_cram(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(archive.to_cram, *args), config)

@require(joint)
def square_batch_region(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(joint.square_batch_region, *args), config)

@require(rnaseq)
def cufflinks_assemble(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(rnaseq.cufflinks_assemble, *args), config)

@require(rnaseq)
def cufflinks_merge(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(rnaseq.cufflinks_merge, *args), config)
