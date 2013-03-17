"""Next-gen variant detection and evaluation with GATK and SnpEff.
"""
import os
import json

from bcbio.log import logger
from bcbio.pipeline.shared import configured_vrn_files
from bcbio.structural import hydra
from bcbio.variation.genotype import variant_filtration, gatk_evaluate_variants
from bcbio.variation import annotation, effects

# ## Genotyping

def postprocess_variants(data):
    """Provide post-processing of variant calls.
    """
    if data["config"]["algorithm"]["snpcall"]:
        logger.info("Finalizing variant calls: %s" % str(data["name"]))
        data["vrn_file"] = finalize_genotyper(data["vrn_file"], data["work_bam"],
                                              data["sam_ref"], data["config"])
        logger.info("Calculating variation effects for %s" % str(data["name"]))
        ann_vrn_file = effects.snpeff_effects(data["vrn_file"], data["genome_build"],
                                              data["config"])
        if ann_vrn_file:
            data["vrn_file"] = ann_vrn_file
    return [[data]]

def finalize_genotyper(call_file, bam_file, ref_file, config):
    """Perform SNP genotyping and analysis.
    """
    vrn_files = configured_vrn_files(config, ref_file)
    variantcaller = config["algorithm"].get("variantcaller", "gatk")
    if variantcaller in ["freebayes", "cortex", "samtools", "gatk-haplotype", "varscan"]:
        call_file = annotation.annotate_nongatk_vcf(call_file, bam_file, vrn_files.dbsnp,
                                                    ref_file, config)
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
