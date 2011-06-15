#!/usr/bin/env python
"""Perform an automated analysis on a sequencing run using Galaxy information.

Given a directory of solexa output, this retrieves details about the sequencing
run from the Galaxy description, and uses this to perform an initial alignment
and analysis.

Usage:
    automated_initial_analysis.py <YAML config file> <flow cell dir>
                                  [<YAML run information>]

The optional <YAML run information> file specifies details about the
flowcell lanes, instead of retrieving it from Galaxy. An example
configuration file is located in 'config/run_info.yaml'

Workflow:
    - Retrieve details on a run.
    - Align fastq files to reference genome.
    - Perform secondary analyses like SNP calling.
    - Generate summary report.
"""
import os
import sys
import json
import subprocess
import glob
import copy
from optparse import OptionParser
import xml.etree.ElementTree as ET
import StringIO

import yaml
import logbook

from bcbio.solexa.flowcell import (get_flowcell_info, get_fastq_dir)
from bcbio.galaxy.api import GalaxyApiAccess
from bcbio import utils
from bcbio.log import create_log_handler
from bcbio.pipeline.fastq import get_fastq_files
from bcbio.pipeline.alignment import align_to_sort_bam, get_genome_ref
from bcbio.pipeline.demultiplex import split_by_barcode, add_multiplex_across_lanes
from bcbio.pipeline.merge import (combine_fastq_files, organize_samples,
                                  merge_bam_files)
from bcbio.pipeline.qcsummary import (generate_align_summary, write_metrics)

LOG_NAME = os.path.splitext(os.path.basename(__file__))[0]
log = logbook.Logger(LOG_NAME)

def main(config_file, fc_dir, run_info_yaml=None):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    log_handler = create_log_handler(config, LOG_NAME)
    
    with log_handler.applicationbound():
        run_main(config, config_file, fc_dir, run_info_yaml)

def run_main(config, config_file, fc_dir, run_info_yaml):
    work_dir = os.getcwd()
    fc_name, fc_date = get_flowcell_info(fc_dir)

    if run_info_yaml and os.path.exists(run_info_yaml):
        log.info("Found YAML samplesheet, using %s instead of Galaxy API" % run_info_yaml)
        with open(run_info_yaml) as in_handle:
            run_details = yaml.load(in_handle)
        run_info = dict(details=run_details, run_id="")
    else:
        log.info("Fetching run details from Galaxy instance")
        galaxy_api = GalaxyApiAccess(config['galaxy_url'], config['galaxy_api_key'])
        run_info = galaxy_api.run_details(fc_name, fc_date)
    fastq_dir = get_fastq_dir(fc_dir)
    run_items = add_multiplex_across_lanes(run_info["details"], fastq_dir, fc_name)
    align_dir = os.path.join(work_dir, "alignments")

    # process each flowcell lane
    with utils.cpmap(config["algorithm"]["num_cores"]) as cpmap:
        for _ in cpmap(process_lane,
                       ((i, fastq_dir, fc_name, fc_date, align_dir, config, config_file)
                        for i in run_items)):
            pass
    # process samples, potentially multiplexed across multiple lanes
    sample_files, sample_fastq, sample_info = organize_samples(align_dir,
            fastq_dir, work_dir, fc_name, fc_date, run_items)
    with utils.cpmap(config["algorithm"]["num_cores"]) as cpmap:
        for _ in cpmap(process_sample, ((name, sample_fastq[name], sample_info[name],
                                         bam_files, work_dir, config, config_file)
                                        for name, bam_files in sample_files)):
            pass
    write_metrics(run_info, work_dir, fc_dir, fc_name, fc_date, fastq_dir)

@utils.map_wrap
def process_lane(info, fastq_dir, fc_name, fc_date, align_dir, config,
        config_file):
    """Do alignments for a lane, potentially splitting based on barcodes.
    """
    config = _update_config_w_custom(config, info)
    
    sample_name = info.get("description", "")
    if (config["algorithm"].get("include_short_name", True) and
            info.get("name", "")):
        sample_name = "%s---%s" % (info.get("name", ""), sample_name)
    genome_build = info.get("genome_build", None)
    multiplex = info.get("multiplex", None)
    
    log.info("Processing sample %s on lane %s with reference genome %s by researcher %s. Using %s analysis preset" \
             % (sample_name, info["lane"], genome_build, \
                info.get("researcher", "unknown"), info.get("analysis", "")))
    if multiplex:
        log.debug("Sample %s is multiplexed as: %s" % (sample_name, multiplex))
    
    fastq_dir, galaxy_dir = _get_full_paths(fastq_dir, config, config_file)
    full_fastq1, full_fastq2 = get_fastq_files(fastq_dir, info['lane'], fc_name)
    lane_name = "%s_%s_%s" % (info['lane'], fc_date, fc_name)
    for mname, msample, fastq1, fastq2 in split_by_barcode(full_fastq1,
            full_fastq2, multiplex, lane_name, config):
        mlane_name = "%s_%s" % (lane_name, mname) if mname else lane_name
        if msample is None:
            msample = "%s---%s" % (sample_name, mname)
        aligner = config["algorithm"].get("aligner", None)
        if os.path.exists(fastq1) and aligner:
            log.info("Aligning lane %s with %s aligner" % (lane_name, aligner))
            align_to_sort_bam(fastq1, fastq2, genome_build, aligner,
                              mlane_name, msample, align_dir, galaxy_dir,
                              config)

@utils.map_wrap
def process_sample(sample_name, fastq_files, info, bam_files, work_dir,
        config, config_file):
    """Finalize processing for a sample, potentially multiplexed.
    """
    config = _update_config_w_custom(config, info)
    
    genome_build = info.get("genome_build", None)
    (_, galaxy_dir) = _get_full_paths("", config, config_file)
    (_, sam_ref) = get_genome_ref(genome_build, config["algorithm"]["aligner"],
                                  galaxy_dir)
    fastq1, fastq2 = combine_fastq_files(fastq_files, work_dir)
    log.info("Combining and preparing wig file %s" % str(sample_name))
    sort_bam = merge_bam_files(bam_files, work_dir, config)
    bam_to_wig(sort_bam, config, config_file)
    if config["algorithm"]["recalibrate"]:
        log.info("Recalibrating %s with GATK" % str(sample_name))
        dbsnp_file = get_dbsnp_file(config, sam_ref)
        gatk_bam = recalibrate_quality(sort_bam, sam_ref,
                dbsnp_file, config_file)
        log.info("Analyzing recalibration %s" % str(sample_name))
        analyze_recalibration(gatk_bam, fastq1, fastq2)
        if config["algorithm"]["snpcall"]:
            log.info("SNP genotyping %s with GATK" % str(sample_name))
            vrn_file = run_genotyper(gatk_bam, sam_ref, dbsnp_file, config_file)
            eval_genotyper(vrn_file, sam_ref, dbsnp_file, config)
            log.info("Calculating variation effects for %s" % str(sample_name))
            variation_effects(vrn_file, genome_build, sam_ref, config)
    if sam_ref is not None:
        print sample_name, "Generating summary files"
        generate_align_summary(sort_bam, fastq2 is not None, sam_ref,
                config, sample_name, config_file)

def bam_to_wig(bam_file, config, config_file):
    """Provide a BigWig coverage file of the sorted alignments.
    """
    wig_file = "%s.bigwig" % os.path.splitext(bam_file)[0]
    if not (os.path.exists(wig_file) and os.path.getsize(wig_file) > 0):
        cl = [config["analysis"]["towig_script"], bam_file, config_file]
        subprocess.check_call(cl)
    return wig_file

def recalibrate_quality(sort_bam_file, sam_ref, dbsnp_file, config_file):
    """Recalibrate alignments with GATK and provide pdf summary.
    """
    bam_file = sort_bam_file.replace("-sort.bam", ".bam")
    cl = ["picard_gatk_recalibrate.py", config_file, sam_ref, bam_file]
    if dbsnp_file:
        cl.append(dbsnp_file)
    subprocess.check_call(cl)
    out_files = glob.glob("%s*-gatkrecal.bam" % os.path.splitext(sort_bam_file)[0])
    assert len(out_files) == 1, out_files
    return out_files[0]

def run_genotyper(bam_file, ref_file, dbsnp_file, config_file):
    """Perform SNP genotyping and analysis using GATK.
    """
    cl = ["gatk_genotyper.py", config_file, ref_file, bam_file]
    if dbsnp_file:
        cl.append(dbsnp_file)
    subprocess.check_call(cl)
    base, ext = os.path.splitext(bam_file)
    return glob.glob("%s*snp-filter.vcf" % base)[0]

def eval_genotyper(vrn_file, ref_file, dbsnp_file, config):
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

def get_dbsnp_file(config, sam_ref):
    snp_file = config["algorithm"].get("dbsnp", None)
    if snp_file:
        base_dir = os.path.dirname(os.path.dirname(sam_ref))
        snp_file = os.path.join(base_dir, snp_file)
    return snp_file

def analyze_recalibration(recal_file, fastq1, fastq2):
    """Provide a pdf report of GATK recalibration of scores.
    """
    cl = ["analyze_quality_recal.py", recal_file, fastq1]
    if fastq2:
        cl.append(fastq2)
    subprocess.check_call(cl)

def _get_full_paths(fastq_dir, config, config_file):
    """Retrieve full paths for directories in the case of relative locations.
    """
    fastq_dir = utils.add_full_path(fastq_dir)
    config_dir = utils.add_full_path(os.path.dirname(config_file))
    galaxy_config_file = utils.add_full_path(config["galaxy_config"], config_dir)
    return fastq_dir, os.path.dirname(galaxy_config_file)

# Utility functions

def _update_config_w_custom(config, lane_info):
    """Update the configuration for this lane if a custom analysis is specified.
    """
    config = copy.deepcopy(config)
    analysis_type = lane_info.get("analysis", "")
    custom = config["custom_algorithms"].get(analysis_type, None)
    if custom:
        for key, val in custom.iteritems():
            config["algorithm"][key] = val
    return config

if __name__ == "__main__":
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    kwargs = dict()
    main(*args, **kwargs)
