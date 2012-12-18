"""Next-gen variant detection and evaluation with GATK and SnpEff.
"""
import os
import json
import subprocess

from bcbio.variation.recalibrate import gatk_recalibrate
from bcbio.variation.genotype import variant_filtration, gatk_evaluate_variants
from bcbio.variation.effects import snpeff_effects
from bcbio.variation.annotation import annotate_effects
from bcbio.variation import freebayes, phasing
from bcbio.pipeline.shared import (configured_vrn_files, configured_ref_file)
from bcbio.structural import hydra

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

def finalize_genotyper(call_file, bam_file, ref_file, config):
    """Perform SNP genotyping and analysis.
    """
    vrn_files = configured_vrn_files(config, ref_file)
    variantcaller = config["algorithm"].get("variantcaller", "gatk")
    if variantcaller in ["freebayes", "cortex", "samtools"]:
        call_file = freebayes.postcall_annotate(call_file, bam_file, ref_file, vrn_files, config)
    filter_snp = variant_filtration(call_file, ref_file, vrn_files, config)
    phase_snp = phasing.read_backed_phasing(filter_snp, bam_file, ref_file, config)
    _eval_genotyper(phase_snp, ref_file, vrn_files.dbsnp, config)
    return phase_snp

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

    Runs snpEff, returning the resulting effects file. No longer runs the GATK
    annotator, since it requires an old version of snpEff.
    """
    return snpeff_effects(vrn_file, genome_build, config)

# ## Structural variation

def detect_sv(data):
    """Detect structural variation for input sample.
    """
    sv_todo = data["config"]["algorithm"].get("sv_detection", None)
    if sv_todo is not None and data.get("fastq2"):
        if sv_todo == "hydra":
            sv_calls = hydra.detect_sv(data["work_bam"], data["genome_build"],
                                       data["dirs"], data["config"])
        else:
            raise ValueError("Unexpected structural variation method:{}".format(sv_todo))
    return [[data]]
