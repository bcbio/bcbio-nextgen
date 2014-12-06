"""Multiprocessing ready entry points for sample analysis.
"""
from bcbio import structural, utils, chipseq
from bcbio.bam import callable
from bcbio.rnaseq import (sailfish, express)
from bcbio.ngsalign import alignprep
from bcbio.pipeline import (archive, disambiguate, qcsummary, sample,
                            main, shared, variation, run_info, rnaseq)
from bcbio.variation import (bamprep, bedutils, coverage, genotype, ensemble,
                             joint, multi, population, recalibrate, validate,
                             vcfutils)

@utils.map_wrap
def run_sailfish(*args):
    return sailfish.run_sailfish(*args)

@utils.map_wrap
def prepare_sample(*args):
    return sample.prepare_sample(*args)

@utils.map_wrap
def trim_sample(*args):
    return sample.trim_sample(*args)

@utils.map_wrap
def process_alignment(*args):
    return sample.process_alignment(*args)

@utils.map_wrap
def postprocess_alignment(*args):
    return sample.postprocess_alignment(*args)

@utils.map_wrap
def prep_samples(*args):
    return sample.prep_samples(*args)

@utils.map_wrap
def prep_align_inputs(*args):
    return alignprep.create_inputs(*args)

@utils.map_wrap
def merge_sample(*args):
    return sample.merge_sample(*args)

@utils.map_wrap
def delayed_bam_merge(*args):
    return sample.delayed_bam_merge(*args)

@utils.map_wrap
def piped_bamprep(*args):
    return bamprep.piped_bamprep(*args)

@utils.map_wrap
def prep_recal(*args):
    return recalibrate.prep_recal(*args)

@utils.map_wrap
def split_variants_by_sample(*args):
    return multi.split_variants_by_sample(*args)

@utils.map_wrap
def postprocess_variants(*args):
    return variation.postprocess_variants(*args)

@utils.map_wrap
def pipeline_summary(*args):
    return qcsummary.pipeline_summary(*args)

@utils.map_wrap
def qsignature_summary(*args):
    return qcsummary.qsignature_summary(*args)

@utils.map_wrap
def generate_transcript_counts(*args):
    return rnaseq.generate_transcript_counts(*args)

@utils.map_wrap
def run_cufflinks(*args):
    return rnaseq.run_cufflinks(*args)

@utils.map_wrap
def run_express(*args):
    return rnaseq.run_express(*args)

@utils.map_wrap
def combine_bam(*args):
    return shared.combine_bam(*args)

@utils.map_wrap
def variantcall_sample(*args):
    return genotype.variantcall_sample(*args)

@utils.map_wrap
def combine_variant_files(*args):
    return vcfutils.combine_variant_files(*args)

@utils.map_wrap
def concat_variant_files(*args):
    return vcfutils.concat_variant_files(*args)

@utils.map_wrap
def merge_variant_files(*args):
    return vcfutils.merge_variant_files(*args)

@utils.map_wrap
def detect_sv(*args):
    return structural.detect_sv(*args)

@utils.map_wrap
def combine_calls(*args):
    return ensemble.combine_calls(*args)

@utils.map_wrap
def prep_gemini_db(*args):
    return population.prep_gemini_db(*args)

@utils.map_wrap
def combine_bed(*args):
    return bedutils.combine(*args)

@utils.map_wrap
def calc_callable_loci(*args):
    return callable.calc_callable_loci(*args)

@utils.map_wrap
def combine_sample_regions(*args):
    return callable.combine_sample_regions(*args)

@utils.map_wrap
def compare_to_rm(*args):
    return validate.compare_to_rm(*args)

@utils.map_wrap
def coverage_summary(*args):
    return coverage.summary(*args)

@utils.map_wrap
def run_disambiguate(*args):
    return disambiguate.run(*args)

@utils.map_wrap
def disambiguate_split(*args):
    return disambiguate.split(*args)

@utils.map_wrap
def clean_chipseq_alignment(*args):
    return chipseq.clean_chipseq_alignment(*args)

@utils.map_wrap
def archive_to_cram(*args):
    return archive.to_cram(*args)

@utils.map_wrap
def square_batch_region(*args):
    return joint.square_batch_region(*args)

@utils.map_wrap
def cufflinks_assemble(*args):
    return rnaseq.cufflinks_assemble(*args)

@utils.map_wrap
def cufflinks_merge(*args):
    return rnaseq.cufflinks_merge(*args)

@utils.map_wrap
def organize_samples(*args):
    return run_info.organize(*args)

@utils.map_wrap
def run_main(*args):
    work_dir, ready_config_file, systemconfig, fcdir, parallel, samples = args
    return main.run_main(work_dir, run_info_yaml=ready_config_file,
                         config_file=systemconfig, fc_dir=fcdir,
                         parallel=parallel, samples=samples)
