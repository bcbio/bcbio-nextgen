"""Next-gen variant detection and evaluation with GATK and SnpEff.
"""
import os
import json
import subprocess

from bcbio.variation.recalibrate import gatk_recalibrate
from bcbio.variation.realign import gatk_realigner
from bcbio.variation.genotype import gatk_genotyper, gatk_evaluate_variants

# ## Recalibration

def recalibrate_quality(sort_bam_file, fastq1, fastq2, sam_ref, config):
    """Recalibrate alignments with GATK and provide pdf summary.
    """
    dbsnp_file = _get_dbsnp_file(config, sam_ref)
    recal_file = gatk_recalibrate(sort_bam_file, sam_ref, config, dbsnp_file)
    _analyze_recalibration(recal_file, fastq1, fastq2)
    return recal_file

def _analyze_recalibration(recal_file, fastq1, fastq2):
    """Provide a pdf report of GATK recalibration of scores.
    """
    cl = ["analyze_quality_recal.py", recal_file, fastq1]
    if fastq2:
        cl.append(fastq2)
    subprocess.check_call(cl)

def _get_dbsnp_file(config, sam_ref):
    snp_file = config["algorithm"].get("dbsnp", None)
    if snp_file:
        base_dir = os.path.dirname(os.path.dirname(sam_ref))
        snp_file = os.path.join(base_dir, snp_file)
    return snp_file

# ## Genotyping

def run_genotyper(bam_file, ref_file, config):
    """Perform SNP genotyping and analysis using GATK.
    """
    dbsnp_file = _get_dbsnp_file(config, ref_file)
    realign_bam = gatk_realigner(bam_file, ref_file, config, dbsnp_file)
    filter_snp = gatk_genotyper(realign_bam, ref_file, config, dbsnp_file)
    _eval_genotyper(filter_snp, ref_file, dbsnp_file, config)
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

def variation_effects(vrn_file, genome_build, ref_file, config):
    """Calculate effects of variations, associating them with transcripts.
    """
    snp_eff_dir = config["program"]["snpEff"]
    snp_eff_jar = os.path.join(snp_eff_dir, "snpEff.jar")
    cl = ["variant_effects.py", snp_eff_jar, vrn_file, genome_build]
    target = config["algorithm"].get("hybrid_target", "")
    if target:
        base_dir = os.path.dirname(os.path.dirname(ref_file))
        cl.append(os.path.join(base_dir, target))
    subprocess.check_call(cl)


