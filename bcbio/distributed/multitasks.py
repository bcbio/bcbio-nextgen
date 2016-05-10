"""Multiprocessing ready entry points for sample analysis.
"""
from bcbio import heterogeneity, hla, structural, utils, chipseq, upload
from bcbio.bam import callable
from bcbio.srna import sample as srna
from bcbio.srna import group as seqcluster
from bcbio.chipseq import peaks
from bcbio.cwl import create as cwl_create
from bcbio.rnaseq import (sailfish, rapmap, salmon, umi)
from bcbio.ngsalign import alignprep
from bcbio.pipeline import (archive, disambiguate, qcsummary, region, sample,
                            main, shared, variation, run_info, rnaseq)
from bcbio.variation import (bamprep, bedutils, genotype, ensemble,
                             joint, multi, population, recalibrate, validate,
                             vcfutils)

@utils.map_wrap
def run_tagcount(*args):
    return umi.tagcount(*args)

@utils.map_wrap
def run_filter_barcodes(*args):
    return umi.filter_barcodes(*args)

@utils.map_wrap
def run_barcode_histogram(*args):
    return umi.barcode_histogram(*args)

@utils.map_wrap
def run_umi_transform(*args):
    return umi.umi_transform(*args)

@utils.map_wrap
def run_salmon_reads(*args):
    return salmon.run_salmon_reads(*args)

@utils.map_wrap
def run_salmon_bam(*args):
    return salmon.run_salmon_bam(*args)

@utils.map_wrap
def run_sailfish(*args):
    return sailfish.run_sailfish(*args)

@utils.map_wrap
def run_rapmap_align(*args):
    return rapmap.run_rapmap_align(*args)

@utils.map_wrap
def prepare_sample(*args):
    return sample.prepare_sample(*args)

@utils.map_wrap
def prepare_bcbio_samples(*args):
    return sample.prepare_bcbio_samples(*args)

@utils.map_wrap
def trim_sample(*args):
    return sample.trim_sample(*args)

@utils.map_wrap
def trim_srna_sample(*args):
    return srna.trim_srna_sample(*args)

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
def srna_annotation(*args):
    return srna.sample_annotation(*args)

@utils.map_wrap
def seqcluster_prepare(*args):
    return seqcluster.run_prepare(*args)

@utils.map_wrap
def seqcluster_cluster(*args):
    return seqcluster.run_cluster(*args)

@utils.map_wrap
def srna_alignment(*args):
    return seqcluster.run_align(*args)

@utils.map_wrap
def peakcalling(*args):
    return peaks.calling(*args)

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
def merge_split_alignments(*args):
    return sample.merge_split_alignments(*args)

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
def multiqc_summary(*args):
    return qcsummary.multiqc_summary(*args)

@utils.map_wrap
def generate_transcript_counts(*args):
    return rnaseq.generate_transcript_counts(*args)

@utils.map_wrap
def run_cufflinks(*args):
    return rnaseq.run_cufflinks(*args)

@utils.map_wrap
def run_stringtie_expression(*args):
    return rnaseq.run_stringtie_expression(*args)

@utils.map_wrap
def run_express(*args):
    return rnaseq.run_express(*args)

@utils.map_wrap
def run_dexseq(*args):
    return rnaseq.run_dexseq(*args)

@utils.map_wrap
def run_rnaseq_variant_calling(*args):
    return rnaseq.run_rnaseq_variant_calling(*args)

@utils.map_wrap
def run_rnaseq_joint_genotyping(*args):
    return rnaseq.run_rnaseq_joint_genotyping(*args)

@utils.map_wrap
def combine_bam(*args):
    return shared.combine_bam(*args)

@utils.map_wrap
def batch_for_variantcall(*args):
    return genotype.batch_for_variantcall(*args)

@utils.map_wrap
def variantcall_batch_region(*args):
    return genotype.variantcall_batch_region(*args)

@utils.map_wrap
def concat_batch_variantcalls(*args):
    return genotype.concat_batch_variantcalls(*args)

@utils.map_wrap
def get_parallel_regions(*args):
    return region.get_parallel_regions(*args)

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
def call_hla(*args):
    return hla.call_hla(*args)

@utils.map_wrap
def detect_sv(*args):
    return structural.detect_sv(*args)

@utils.map_wrap
def validate_sv(*args):
    return structural.validate_sv(*args)

@utils.map_wrap
def heterogeneity_estimate(*args):
    return heterogeneity.estimate(*args)

@utils.map_wrap
def finalize_sv(*args):
    return structural.finalize_sv(*args)

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
def run_disambiguate(*args):
    return disambiguate.run(*args)

@utils.map_wrap
def disambiguate_split(*args):
    return disambiguate.split(*args)

@utils.map_wrap
def disambiguate_merge_extras(*args):
    return disambiguate.merge_extras(*args)

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
def stringtie_merge(*args):
    return rnaseq.stringtie_merge(*args)

@utils.map_wrap
def organize_samples(*args):
    return run_info.organize(*args)

@utils.map_wrap
def prep_system(*args):
    return run_info.prep_system(*args)

@utils.map_wrap
def upload_samples(*args):
    return upload.from_sample(*args)

@utils.map_wrap
def upload_samples_project(*args):
    return upload.project_from_sample(*args)

@utils.map_wrap
def create_cwl(*args):
    return cwl_create.from_world(*args)

@utils.map_wrap
def run_main(*args):
    work_dir, ready_config_file, systemconfig, fcdir, parallel, samples = args
    return main.run_main(work_dir, run_info_yaml=ready_config_file,
                         config_file=systemconfig, fc_dir=fcdir,
                         parallel=parallel, samples=samples)
