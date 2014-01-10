"""Multiprocessing ready entry points for sample analysis.
"""
from bcbio import structural, utils, chipseq
from bcbio.bam import callable
from bcbio.ngsalign import alignprep, tophat, star
from bcbio.pipeline import (disambiguate, lane, qcsummary, sample, shared, variation,
                            rnaseq, alignment)
from bcbio.variation import (bamprep, coverage, realign, genotype, ensemble, multi, population,
                             recalibrate, validate, vcfutils)

@utils.map_wrap
def process_lane(*args):
    return lane.process_lane(*args)

@utils.map_wrap
def trim_lane(*args):
    return lane.trim_lane(*args)

@utils.map_wrap
def process_alignment(*args):
    return lane.process_alignment(*args)
process_alignment.metadata = {"resources": ["novoalign", "bwa", "bowtie2", "tophat2"],
                              "ensure": {"tophat2": tophat.job_requirements}}

@utils.map_wrap
def postprocess_alignment(*args):
    return lane.postprocess_alignment(*args)
postprocess_alignment.metadata = {"resources": ["gatk"]}

@utils.map_wrap
def prep_align_inputs(*args):
    return alignprep.create_inputs(*args)

@utils.map_wrap
def merge_sample(*args):
    return sample.merge_sample(*args)

@utils.map_wrap
def delayed_bam_merge(*args):
    return sample.delayed_bam_merge(*args)
delayed_bam_merge.metadata = {"resources": ["samtools"]}

@utils.map_wrap
def piped_bamprep(*args):
    return bamprep.piped_bamprep(*args)

@utils.map_wrap
def prep_recal(*args):
    return recalibrate.prep_recal(*args)
prep_recal.metadata = {"resources": ["gatk"]}

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
    return variation.postprocess_variants(*args)
postprocess_variants.metadata = {"resources": ["gatk-vqsr", "gatk", "snpEff"]}

@utils.map_wrap
def pipeline_summary(*args):
    return qcsummary.pipeline_summary(*args)

@utils.map_wrap
def generate_transcript_counts(*args):
    return rnaseq.generate_transcript_counts(*args)

@utils.map_wrap
def run_cufflinks(*args):
    return rnaseq.run_cufflinks(*args)
run_cufflinks.metadata = {"resources": ["cufflinks"]}

@utils.map_wrap
def generate_bigwig(*args):
    return sample.generate_bigwig(*args)

@utils.map_wrap
def combine_bam(*args):
    return shared.combine_bam(*args)

@utils.map_wrap
def variantcall_sample(*args):
    return genotype.variantcall_sample(*args)
variantcall_sample.metadata = {"resources": ["gatk", "freebayes", "gatk-haplotype"]}

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
prep_gemini_db.metadata = {"resources": ["gemini"]}

@utils.map_wrap
def combine_bed(*args):
    return callable.combine_bed(*args)

@utils.map_wrap
def calc_callable_loci(*args):
    return callable.calc_callable_loci(*args)

@utils.map_wrap
def compare_to_rm(*args):
    return validate.compare_to_rm(*args)

@utils.map_wrap
def coverage_summary(*args):
    return coverage.summary(*args)
coverage_summary.metadata = {"resources": ["bcbio_coverage"]}

@utils.map_wrap
def run_disambiguate(*args):
    return disambiguate.run(*args)

@utils.map_wrap
def clean_chipseq_alignment(*args):
    return chipseq.clean_chipseq_alignment(*args)
