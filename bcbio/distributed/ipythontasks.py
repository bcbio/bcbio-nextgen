"""Ipython parallel ready entry points for parallel execution
"""
import contextlib
import os

try:
    from ipyparallel import require
except ImportError:
    from IPython.parallel import require

from bcbio import heterogeneity, hla, chipseq, structural, upload
from bcbio.bam import callable
from bcbio.rnaseq import (sailfish, rapmap, salmon, umi)
from bcbio.distributed import ipython
from bcbio.ngsalign import alignprep
from bcbio import rnaseq
from bcbio.srna import sample as srna
from bcbio.srna import group as seqcluster
from bcbio.chipseq import peaks
from bcbio.pipeline import (archive, config_utils, disambiguate, sample,
                            qcsummary, shared, variation, run_info, rnaseq)
from bcbio.provenance import system
from bcbio.variation import (bamprep, coverage, genotype, ensemble, joint,
                             multi, population, recalibrate, validate, vcfutils)
from bcbio.log import logger, setup_local_logging

@contextlib.contextmanager
def _setup_logging(args):
    # Set environment to standard to use periods for decimals and avoid localization
    os.environ["LC_ALL"] = "C"
    os.environ["LC"] = "C"
    os.environ["LANG"] = "C"
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

# Potential wrapper to avoid boilerplate if we can get dill working for closures
from functools import wraps
def _pack_n_log(f):
    from bcbio.distributed import ipython
    @wraps(f)
    def wrapper(*args):
        args = ipython.unzip_args(args)
        with _setup_logging(args) as config:
            return ipython.zip_args(fn(*args))
    return wrapper

@require(sample)
def prepare_sample(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(sample.prepare_sample, *args))

@require(sample)
def prepare_bcbio_samples(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(sample.prepare_bcbio_samples, *args))

@require(sample)
def trim_sample(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(sample.trim_sample, *args))

@require(srna)
def trim_srna_sample(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(srna.trim_srna_sample, *args))

@require(srna)
def srna_annotation(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(srna.sample_annotation, *args))

@require(seqcluster)
def seqcluster_prepare(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(seqcluster.run_prepare, *args))

@require(seqcluster)
def seqcluster_cluster(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(seqcluster.run_cluster, *args))

@require(seqcluster)
def srna_alignment(* args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(seqcluster.run_align, *args))

@require(peaks)
def peakcalling(* args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(peaks.calling, *args))

@require(sailfish)
def run_sailfish(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args):
        return ipython.zip_args(apply(sailfish.run_sailfish, *args))

@require(rapmap)
def run_rapmap_align(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args):
        return ipython.zip_args(apply(rapmap.run_rapmap_align, *args))

@require(umi)
def run_umi_transform(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args):
        return ipython.zip_args(apply(umi.umi_transform, *args))

@require(umi)
def run_filter_barcodes(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args):
        return ipython.zip_args(apply(umi.filter_barcodes, *args))

@require(umi)
def run_tagcount(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args):
        return ipython.zip_args(apply(umi.tagcount, *args))

@require(umi)
def run_barcode_histogram(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args):
        return ipython.zip_args(apply(umi.barcode_histogram, *args))

@require(salmon)
def run_salmon_bam(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args):
        return ipython.zip_args(apply(salmon.run_salmon_bam, *args))

@require(salmon)
def run_salmon_reads(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args):
        return ipython.zip_args(apply(salmon.run_salmon_reads, *args))

@require(sample)
def process_alignment(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(sample.process_alignment, *args))

@require(alignprep)
def prep_align_inputs(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(alignprep.create_inputs, *args))

@require(sample)
def postprocess_alignment(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(sample.postprocess_alignment, *args))

@require(sample)
def prep_samples(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(sample.prep_samples, *args))

@require(sample)
def merge_sample(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(sample.merge_sample, *args))

@require(sample)
def delayed_bam_merge(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(sample.delayed_bam_merge, *args))

@require(sample)
def merge_split_alignments(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(sample.merge_split_alignments, *args))

@require(sample)
def recalibrate_sample(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(sample.recalibrate_sample, *args))

@require(recalibrate)
def prep_recal(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(recalibrate.prep_recal, *args))

@require(multi)
def split_variants_by_sample(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(multi.split_variants_by_sample, *args))

@require(bamprep)
def piped_bamprep(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(bamprep.piped_bamprep, *args))

@require(variation)
def postprocess_variants(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(variation.postprocess_variants, *args))

@require(qcsummary)
def pipeline_summary(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(qcsummary.pipeline_summary, *args))

@require(qcsummary)
def qsignature_summary(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(qcsummary.qsignature_summary, *args))

@require(qcsummary)
def multiqc_summary(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(qcsummary.multiqc_summary, *args))

@require(rnaseq)
def generate_transcript_counts(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(rnaseq.generate_transcript_counts, *args))

@require(rnaseq)
def run_cufflinks(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(rnaseq.run_cufflinks, *args))

@require(rnaseq)
def run_stringtie_expression(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(rnaseq.run_stringtie_expression, *args))

@require(rnaseq)
def run_rnaseq_variant_calling(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(rnaseq.run_rnaseq_variant_calling, *args))

@require(rnaseq)
def run_rnaseq_joint_genotyping(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(rnaseq.run_rnaseq_joint_genotyping, *args))

@require(rnaseq)
def run_express(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(rnaseq.run_express, *args))

@require(rnaseq)
def run_dexseq(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(rnaseq.run_dexseq, *args))


@require(shared)
def combine_bam(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(shared.combine_bam, *args))

@require(callable)
def combine_sample_regions(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(callable.combine_sample_regions, *args))

@require(genotype)
def variantcall_sample(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(genotype.variantcall_sample, *args))

@require(vcfutils)
def combine_variant_files(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(vcfutils.combine_variant_files, *args))

@require(vcfutils)
def concat_variant_files(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(vcfutils.concat_variant_files, *args))

@require(vcfutils)
def merge_variant_files(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(vcfutils.merge_variant_files, *args))

@require(population)
def prep_gemini_db(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(population.prep_gemini_db, *args))

@require(hla)
def call_hla(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(hla.call_hla, *args))

@require(structural)
def detect_sv(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(structural.detect_sv, *args))

@require(structural)
def validate_sv(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(structural.validate_sv, *args))

@require(structural)
def finalize_sv(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(structural.finalize_sv, *args))

@require(heterogeneity)
def heterogeneity_estimate(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(heterogeneity.estimate, *args))

@require(ensemble)
def combine_calls(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(ensemble.combine_calls, *args))

@require(validate)
def compare_to_rm(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(validate.compare_to_rm, *args))

@require(disambiguate)
def run_disambiguate(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(disambiguate.run, *args))

@require(disambiguate)
def disambiguate_split(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(disambiguate.split, *args))

@require(disambiguate)
def disambiguate_merge_extras(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(disambiguate.merge_extras, *args))

@require(system)
def machine_info(*args):
    args = ipython.unzip_args(args)
    return ipython.zip_args(system.machine_info())

@require(chipseq)
def clean_chipseq_alignment(*args):
    args = ipython.unzip_args(args)
    return ipython.zip_args(apply(chipseq.clean_chipseq_alignment, *args))

@require(archive)
def archive_to_cram(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(archive.to_cram, *args))

@require(joint)
def square_batch_region(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(joint.square_batch_region, *args))

@require(rnaseq)
def cufflinks_assemble(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(rnaseq.cufflinks_assemble, *args))

@require(rnaseq)
def cufflinks_merge(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(rnaseq.cufflinks_merge, *args))

@require(rnaseq)
def stringtie_merge(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(rnaseq.stringtie_merge, *args))

@require(run_info)
def organize_samples(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(run_info.organize, *args))

@require(run_info)
def prep_system(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(run_info.prep_system, *args))

@require(upload)
def upload_samples(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(upload.from_sample, *args))

@require(upload)
def upload_samples_project(*args):
    args = ipython.unzip_args(args)
    with _setup_logging(args) as config:
        return ipython.zip_args(apply(upload.project_from_sample, *args))
