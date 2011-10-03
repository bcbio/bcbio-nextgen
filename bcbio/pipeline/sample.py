"""High level entry point for processing a sample.

Samples may include multiple lanes, or barcoded subsections of lanes,
processed together.
"""
import os
import subprocess

from bcbio.utils import file_transaction
from bcbio.pipeline.lane import _update_config_w_custom
from bcbio.pipeline import log
from bcbio.pipeline.alignment import get_genome_ref
from bcbio.pipeline.merge import (combine_fastq_files, merge_bam_files)
from bcbio.pipeline.qcsummary import generate_align_summary
from bcbio.pipeline.variation import (recalibrate_quality, run_genotyper,
                                      variation_effects, configured_ref_file)
from bcbio.variation.realign import gatk_realigner
from bcbio.rnaseq.cufflinks import assemble_transcripts

def merge_sample(sample_name, fastq_files, info, bam_files, dirs,
                 config, config_file):
    """Merge fastq and BAM files for multiple samples.
    """
    log.info("Combining fastq and BAM files %s" % str(sample_name))
    config = _update_config_w_custom(config, info)
    fastq1, fastq2 = combine_fastq_files(fastq_files, dirs["work"])
    sort_bam = merge_bam_files(bam_files, dirs["work"], config)
    return [(sample_name, sort_bam, fastq1, fastq2, info, dirs,
             config, config_file)]

def _ref_genome_info(info, config, dirs):
    genome_build = info.get("genome_build", None)
    (_, sam_ref) = get_genome_ref(genome_build, config["algorithm"]["aligner"],
                                  dirs["galaxy"])
    return genome_build, sam_ref

def recalibrate_sample(sample_name, sort_bam, fastq1, fastq2, info,
                       dirs, config, config_file):
    """Recalibrate quality values from aligned sample BAM file.
    """
    log.info("Recalibrating %s with GATK" % str(sample_name))
    _, sam_ref = _ref_genome_info(info, config, dirs)
    if config["algorithm"]["recalibrate"]:
        gatk_bam = recalibrate_quality(sort_bam, fastq1, fastq2, sam_ref,
                                       dirs, config)
    else:
        gatk_bam = sort_bam
    return [(sample_name, gatk_bam, fastq1, fastq2, info,
             dirs, config, config_file)]

def realign_sample(sample_name, bam_file, fastq1, fastq2, info,
                   dirs, config, config_file):
    """Realign sample BAM file at indels.
    """
    log.info("Realigning %s with GATK" % str(sample_name))
    _, sam_ref = _ref_genome_info(info, config, dirs)
    if config["algorithm"]["snpcall"]:
        realign_bam = gatk_realigner(bam_file, sam_ref, config,
                                     configured_ref_file("dbsnp", config, sam_ref))
    else:
        realign_bam = bam_file
    return [(sample_name, realign_bam, fastq1, fastq2, info,
             dirs, config, config_file)]

def process_sample(sample_name, bam_file, fastq1, fastq2, info,
                   dirs, config, config_file):
    """Finalize processing for a sample, potentially multiplexed.
    """
    genome_build = info.get("genome_build", None)
    (_, sam_ref) = get_genome_ref(genome_build, config["algorithm"]["aligner"],
                                  dirs["galaxy"])
    (vrn_file, effects_file) = ("", "")
    if config["algorithm"]["snpcall"]:
        log.info("SNP genotyping %s with GATK" % str(sample_name))
        vrn_file = run_genotyper(bam_file, sam_ref, config)
        log.info("Calculating variation effects for %s" % str(sample_name))
        annotated_vrn_file, effects_file = variation_effects(vrn_file, sam_ref,
                                                             genome_build, config)
    if config["algorithm"].get("transcript_assemble", False):
        tx_file = assemble_transcripts(bam_file, sam_ref, config)
    if sam_ref is not None:
        log.info("Generating summary files: %s" % str(sample_name))
        generate_align_summary(bam_file, fastq2 is not None, sam_ref,
                               sample_name, config, dirs)
    log.info("Preparing wig file %s" % str(sample_name))
    bam_to_wig(bam_file, config, config_file)
    return [(sample_name, fastq1, fastq2, info, bam_file,
             annotated_vrn_file if annotated_vrn_file else vrn_file,
             effects_file)]

def bam_to_wig(bam_file, config, config_file):
    """Provide a BigWig coverage file of the sorted alignments.
    """
    wig_file = "%s.bigwig" % os.path.splitext(bam_file)[0]
    if not (os.path.exists(wig_file) and os.path.getsize(wig_file) > 0):
        cl = [config["analysis"]["towig_script"], bam_file, config_file]
        with file_transaction(wig_file):
            subprocess.check_call(cl)
    return wig_file

