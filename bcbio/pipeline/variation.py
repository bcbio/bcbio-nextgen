"""Next-gen variant detection and evaluation with GATK and SnpEff.
"""
import os
import json
import subprocess

from bcbio.variation.genotype import variant_filtration, gatk_evaluate_variants
from bcbio.variation.effects import snpeff_effects
from bcbio.variation.annotation import annotate_effects
from bcbio.variation import freebayes, phasing
from bcbio.pipeline.shared import (configured_vrn_files, configured_ref_file)
from bcbio.structural import hydra

# ## Genotyping

def finalize_genotyper(call_file, bam_file, ref_file, config):
    """Perform SNP genotyping and analysis.
    """
    vrn_files = configured_vrn_files(config, ref_file)
    variantcaller = config["algorithm"].get("variantcaller", "gatk")
    if variantcaller in ["freebayes", "cortex", "samtools", "gatk-haplotype", "varscan"]:
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
