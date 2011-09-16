"""Next-gen variant detection and evaluation with GATK and SnpEff.
"""
import os
import json
import subprocess

from bcbio.variation.recalibrate import gatk_recalibrate
from bcbio.variation.realign import gatk_realigner
from bcbio.variation.genotype import gatk_genotyper, gatk_evaluate_variants
from bcbio.variation.effects import snpeff_effects

# ## Recalibration

def recalibrate_quality(sort_bam_file, fastq1, fastq2, sam_ref,
                        dirs, config):
    """Recalibrate alignments with GATK and provide pdf summary.
    """
    dbsnp_file = _configured_ref_file("dbsnp", config, ref_file)
    recal_file = gatk_recalibrate(sort_bam_file, sam_ref, config, dbsnp_file)
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

def _configured_ref_file(name, config, sam_ref):
    """Full path to a reference file specified in the configuration.

    Resolves non-absolute paths relative to the base genome reference directory.
    """
    ref_file = config["algorithm"].get(name, None)
    if ref_file:
        if not os.path.isabs(ref_file):
            base_dir = os.path.dirname(os.path.dirname(sam_ref))
            ref_file = os.path.join(base_dir, ref_file)
    return ref_file

# ## Genotyping

def run_genotyper(bam_file, ref_file, config):
    """Perform SNP genotyping and analysis using GATK.
    """
    dbsnp_file = _configured_ref_file("dbsnp", config, ref_file)
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

def variation_effects(vrn_file, genome_build, config):
    """Calculate effects of variations, associating them with transcripts.
    """
    snpeff_jar = os.path.join(config["program"]["snpEff"], "snpEff.jar")
    java_memory = config["algorithm"].get("java_memory", None)
    return snpeff_effects(snpeff_jar, vrn_file, genome_build,
                          config["algorithm"].get("hybrid_target", None),
                          java_memory)
