"""High level entry point for processing a sample.

Samples may include multiple lanes, or barcoded subsections of lanes,
that should be processed together.
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
                                      variation_effects)

def process_sample(sample_name, fastq_files, info, bam_files, dirs,
                   config, config_file):
    """Finalize processing for a sample, potentially multiplexed.
    """
    config = _update_config_w_custom(config, info)

    genome_build = info.get("genome_build", None)
    (_, sam_ref) = get_genome_ref(genome_build, config["algorithm"]["aligner"],
                                  dirs["galaxy"])
    fastq1, fastq2 = combine_fastq_files(fastq_files, dirs["work"])
    log.info("Combining and preparing wig file %s" % str(sample_name))
    sort_bam = merge_bam_files(bam_files, dirs["work"], config)
    bam_to_wig(sort_bam, config, config_file)
    (gatk_bam, vrn_file, effects_file) = ("", "", "")
    if config["algorithm"]["recalibrate"]:
        log.info("Recalibrating %s with GATK" % str(sample_name))
        gatk_bam = recalibrate_quality(sort_bam, fastq1, fastq2, sam_ref, config)
        if config["algorithm"]["snpcall"]:
            log.info("SNP genotyping %s with GATK" % str(sample_name))
            vrn_file = run_genotyper(gatk_bam, sam_ref, config)
            log.info("Calculating variation effects for %s" % str(sample_name))
            effects_file = variation_effects(vrn_file, genome_build, config)
    if sam_ref is not None:
        log.info("Generating summary files: %s" % str(sample_name))
        generate_align_summary(sort_bam, fastq2 is not None, sam_ref,
                config, sample_name, config_file)
    return [sample_name, fastq_files, info, sort_bam, gatk_bam, vrn_file,
            effects_file]

def bam_to_wig(bam_file, config, config_file):
    """Provide a BigWig coverage file of the sorted alignments.
    """
    wig_file = "%s.bigwig" % os.path.splitext(bam_file)[0]
    if not (os.path.exists(wig_file) and os.path.getsize(wig_file) > 0):
        cl = [config["analysis"]["towig_script"], bam_file, config_file]
        with file_transaction(wig_file):
            subprocess.check_call(cl)
    return wig_file

