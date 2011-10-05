"""Next-gen variant detection and evaluation with GATK and SnpEff.
"""
import os
import json
import subprocess

from bcbio.variation.recalibrate import gatk_recalibrate
from bcbio.variation.genotype import variant_filtration, gatk_evaluate_variants
from bcbio.variation.effects import snpeff_effects
from bcbio.variation.annotation import annotate_effects
from bcbio.pipeline.shared import (configured_vrn_files, configured_ref_file)

# ## Recalibration

def recalibrate_quality(sort_bam_file, fastq1, fastq2, sam_ref,
                        dirs, config):
    """Recalibrate alignments with GATK and provide pdf summary.
    """
    dbsnp_file = configured_ref_file("dbsnp", config, sam_ref)
    recal_file = gatk_recalibrate(sort_bam_file, sam_ref, config, dbsnp_file)
    if config["algorithm"].get("recalibration_plots", False):
        _analyze_recalibration(recal_file, fastq1, fastq2, dirs, config)
    return recal_file

def _analyze_recalibration(recal_file, fastq1, fastq2, dirs, config):
    """Provide a pdf report of GATK recalibration of scores.
    """
    qual_opts = {"illumina": "fastq-illumina", "standard": "fastq"}
    qual_format = config["algorithm"].get("quality_format", "illumina").lower()
    cl = ["analyze_quality_recal.py", recal_file, fastq1]
    if fastq2:
        cl.append(fastq2)
    cl.append("--workdir=%s" % dirs["work"])
    cl.append("--input_format=%s" % qual_opts[qual_format])
    subprocess.check_call(cl)

# ## Genotyping

def finalize_genotyper(call_file, ref_file, config):
    """Perform SNP genotyping and analysis using GATK.
    """
    vrn_files = configured_vrn_files(config, ref_file)
    filter_snp = variant_filtration(call_file, ref_file, vrn_files, config)
    _eval_genotyper(filter_snp, ref_file, vrn_files.dbsnp, config)
    return filter_snp

def _eval_genotyper(vrn_file, ref_file, dbsnp_file, config):
    """Evaluate variant genotyping, producing a JSON metrics file with values.
    """
    metrics_file = "%s.eval_metrics" % vrn_file
    target = config["algorithm"].get("hybrid_target", None)
    if not os.path.exists(metrics_file):
        stats = gatk_evaluate_variants(vrn_file, ref_file, config, dbsnp_file, target)
        with open(metrics_file, "w") as out_handle:
            json.dump(stats, out_handle)
    return metrics_file

# ## Calculate variation effects

def variation_effects(vrn_file, genome_file, genome_build, config):
    """Calculate effects of variations, associating them with transcripts.
    """
    snpeff_vcf, snpeff_txt = snpeff_effects(vrn_file, genome_build, config)
    annotated_vcf = annotate_effects(vrn_file, snpeff_vcf, genome_file, config) \
                    if snpeff_vcf else None
    return annotated_vcf, snpeff_txt
