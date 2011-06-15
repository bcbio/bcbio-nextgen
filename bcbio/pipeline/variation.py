"""Next-gen variant detection and evaluation with GATK and SnpEff.
"""
import os
import glob
import subprocess

# ## Recalibration

def recalibrate_quality(sort_bam_file, fastq1, fastq2, sam_ref, config, config_file):
    """Recalibrate alignments with GATK and provide pdf summary.
    """
    dbsnp_file = _get_dbsnp_file(config, sam_ref)
    bam_file = sort_bam_file.replace("-sort.bam", ".bam")
    cl = ["picard_gatk_recalibrate.py", config_file, sam_ref, bam_file]
    if dbsnp_file:
        cl.append(dbsnp_file)
    subprocess.check_call(cl)
    out_files = glob.glob("%s*-gatkrecal.bam" % os.path.splitext(sort_bam_file)[0])
    assert len(out_files) == 1, out_files
    recal_file = out_files[0]
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

def run_genotyper(bam_file, ref_file, config, config_file):
    """Perform SNP genotyping and analysis using GATK.
    """
    dbsnp_file = _get_dbsnp_file(config, ref_file)
    cl = ["gatk_genotyper.py", config_file, ref_file, bam_file]
    if dbsnp_file:
        cl.append(dbsnp_file)
    subprocess.check_call(cl)
    base = os.path.splitext(bam_file)[0]
    vrn_file = glob.glob("%s*snp-filter.vcf" % base)[0]
    _eval_genotyper(vrn_file, ref_file, dbsnp_file, config)
    return vrn_file

def _eval_genotyper(vrn_file, ref_file, dbsnp_file, config):
    """Evaluate variant genotyping, producing a JSON metrics file with values.
    """
    metrics_file = "%s.eval_metrics" % vrn_file
    cl = ["gatk_variant_eval.py", config["program"].get("gatk", config["program"]["picard"]),
          vrn_file, ref_file, dbsnp_file]
    target = config["algorithm"].get("hybrid_target", "")
    if target:
        base_dir = os.path.dirname(os.path.dirname(ref_file))
        cl.append(os.path.join(base_dir, target))
    with open(metrics_file, "w") as out_handle:
        subprocess.check_call(cl, stdout=out_handle)

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


