"""Multiprocessing ready entry points for sample analysis.
"""
from bcbio import heterogeneity, hla, structural, utils, chipseq, upload
from bcbio.bam import callable
from bcbio.srna import sample as srna
from bcbio.srna import group as seqcluster
from bcbio.chipseq import peaks
from bcbio.wgbsseq import cpg_caller, deduplication, trimming
from bcbio.cwl import create as cwl_create
from bcbio.cwl import cwlutils
from bcbio.rnaseq import (sailfish, rapmap, salmon, umi, kallisto, spikein,
                          bcbiornaseq)
from bcbio.ngsalign import alignprep
from bcbio.pipeline import (archive, alignment, disambiguate, qcsummary, region, sample,
                            main, shared, variation, run_info, rnaseq)
from bcbio.qc import multiqc, qsignature
from bcbio.structural import regions as svregions
from bcbio.variation import (bamprep, genotype, ensemble,
                             joint, multi, population, validate,
                             vcfutils, peddy)

@utils.map_wrap
def run_peddy(*args):
    return peddy.run_peddy(*args)

@utils.map_wrap
def run_tagcount(*args):
    return umi.tagcount(*args)

@utils.map_wrap
def run_concatenate_sparse_counts(*args):
    return umi.concatenate_sparse_counts(*args)

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
def demultiplex_samples(*args):
    return umi.demultiplex_samples(*args)

@utils.map_wrap
def run_kallisto_singlecell(*args):
    return kallisto.run_kallisto_singlecell(*args)

@utils.map_wrap
def run_kallisto_index(*args):
    return kallisto.run_kallisto_index(*args)

@utils.map_wrap
def run_kallisto_rnaseq(*args):
    return kallisto.run_kallisto_rnaseq(*args)

@utils.map_wrap
def run_salmon_decoy(*args):
    return salmon.run_salmon_decoy(*args)

@utils.map_wrap
def run_salmon_reads(*args):
    return salmon.run_salmon_reads(*args)

@utils.map_wrap
def run_salmon_bam(*args):
    return salmon.run_salmon_bam(*args)

@utils.map_wrap
def run_salmon_index(*args):
    return salmon.run_salmon_index(*args)

@utils.map_wrap
def run_rapmap_index(*args):
    return rapmap.run_rapmap_index(*args)

@utils.map_wrap
def run_counts_spikein(*args):
    return spikein.run_counts_spikein(*args)

@utils.map_wrap
def run_bcbiornaseqload(*args):
    return bcbiornaseq.make_bcbiornaseq_object(*args)

@utils.map_wrap
def run_sailfish(*args):
    return sailfish.run_sailfish(*args)

@utils.map_wrap
def run_sailfish_index(*args):
    return sailfish.run_sailfish_index(*args)

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
def process_alignment_to_rec(*args):
    return cwlutils.to_rec(*args)

@utils.map_wrap
def process_alignment(*args):
    return sample.process_alignment(*args)

@utils.map_wrap
def alignment_to_rec(*args):
    default_keys = ["config__algorithm__align_split_size",
                    "config__algorithm__aligner",
                    "config__algorithm__mark_duplicates",
                    "reference__bwa__indexes",
                    "reference__snap__indexes",
                    "reference__bowtie2__indexes",
                    "reference__novoalign__indexes",
                    "rgnames__pl", "rgnames__sample", "rgnames__pu",
                    "rgnames__lane", "rgnames__rg", "rgnames__lb"]
    return cwlutils.to_rec_single(*args, default_keys=default_keys)

@utils.map_wrap
def organize_noalign(*args):
    return alignment.organize_noalign(args)

@utils.map_wrap
def postprocess_alignment_to_rec(*args):
    default_keys = ["config__algorithm__coverage_interval", "config__algorithm__seq2c_bed_ready",
                    "config__algorithm__coverage", "config__algorithm__coverage_merged",
                    "config__algorithm__coverage_orig", "config__algorithm__variant_regions",
                    "config__algorithm__variant_regions_merged", "config__algorithm__variant_regions_orig"]
    return cwlutils.to_rec(*args, default_keys=default_keys)

@utils.map_wrap
def postprocess_alignment(*args):
    return sample.postprocess_alignment(*args)

@utils.map_wrap
def prep_samples(*args):
    return sample.prep_samples(*args)

@utils.map_wrap
def prep_samples_to_rec(*args):
    return cwlutils.to_rec(*args)

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
def trim_bs_sample(*args):
    return trimming.trim(*args)


@utils.map_wrap
def cpg_calling(*args):
    return cpg_caller.calling(*args)


@utils.map_wrap
def cpg_processing(*args):
    return cpg_caller.cpg_postprocessing(*args)


@utils.map_wrap
def cpg_stats(*args):
    return cpg_caller.cpg_stats(*args)


@utils.map_wrap
def deduplicate_bismark(*args):
    return deduplication.dedup_bismark(*args)


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
def split_variants_by_sample(*args):
    return multi.split_variants_by_sample(*args)

@utils.map_wrap
def postprocess_variants(*args):
    return variation.postprocess_variants(*args)

@utils.map_wrap
def pipeline_summary(*args):
    return qcsummary.pipeline_summary(*args)

@utils.map_wrap
def qc_to_rec(*args):
    return qcsummary.qc_to_rec(*args)

@utils.map_wrap
def qsignature_summary(*args):
    return qsignature.summary(*args)

@utils.map_wrap
def multiqc_summary(*args):
    return multiqc.summary(*args)

@utils.map_wrap
def generate_transcript_counts(*args):
    return rnaseq.generate_transcript_counts(*args)

@utils.map_wrap
def detect_fusions(*args):
    return rnaseq.detect_fusions(*args)

@utils.map_wrap
def rnaseq_quantitate(*args):
    return rnaseq.quantitate(*args)

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
def run_rnaseq_ann_filter(*args):
    return rnaseq.run_rnaseq_ann_filter(*args)

@utils.map_wrap
def combine_bam(*args):
    return shared.combine_bam(*args)

@utils.map_wrap
def batch_for_variantcall(*args):
    return genotype.batch_for_variantcall(*args)

@utils.map_wrap
def vc_output_record(*args):
    return genotype.vc_output_record(*args)

@utils.map_wrap
def variantcall_batch_region(*args):
    return genotype.variantcall_batch_region(*args)

@utils.map_wrap
def concat_batch_variantcalls(*args):
    return genotype.concat_batch_variantcalls(*args)

@utils.map_wrap
def get_parallel_regions(*args):
    return region.get_parallel_regions_block(*args)

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
def hla_to_rec(*args):
    return cwlutils.to_rec(*args)

@utils.map_wrap
def call_hla(*args):
    return hla.call_hla(*args)

@utils.map_wrap
def calculate_sv_bins(*args):
    return svregions.calculate_sv_bins(*args)

@utils.map_wrap
def calculate_sv_coverage(*args):
    return svregions.calculate_sv_coverage(*args)

@utils.map_wrap
def normalize_sv_coverage(*args):
    return svregions.normalize_sv_coverage(*args)

@utils.map_wrap
def batch_for_sv(*args):
    return structural.batch_for_sv(*args)

@utils.map_wrap
def detect_sv(*args):
    return structural.detect_sv(*args)

@utils.map_wrap
def summarize_sv(*args):
    return structural.summarize_sv(*args)

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
def batch_for_ensemble(*args):
    return ensemble.batch(*args)

@utils.map_wrap
def prep_gemini_db(*args):
    return population.prep_gemini_db(*args)

@utils.map_wrap
def combine_sample_regions(*args):
    return callable.combine_sample_regions(*args)

@utils.map_wrap
def compare_to_rm(*args):
    return validate.compare_to_rm(*args)

@utils.map_wrap
def summarize_vc(*args):
    return variation.summarize_vc(*args)

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
def batch_for_jointvc(*args):
    return joint.batch_for_jointvc(*args)

@utils.map_wrap
def run_jointvc(*args):
    return joint.run_jointvc(*args)

@utils.map_wrap
def finalize_jointvc(*args):
    return joint.finalize_jointvc(*args)

@utils.map_wrap
def get_parallel_regions_jointvc(*args):
    return region.get_parallel_regions(*args)

@utils.map_wrap
def concat_batch_variantcalls_jointvc(*args):
    return joint.concat_batch_variantcalls_jointvc(*args)

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
