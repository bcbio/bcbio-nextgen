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
import subprocess
import copy
from optparse import OptionParser

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
from bcbio.pipeline.variation import (recalibrate_quality, run_genotyper,
                                      variation_effects)

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
    align_dir = os.path.join(work_dir, "alignments")

    fc_name, fc_date = get_flowcell_info(fc_dir)
    run_info = _get_run_info(fc_name, fc_date, config, run_info_yaml)
    fastq_dir, galaxy_dir = _get_full_paths(get_fastq_dir(fc_dir), config, config_file)
    dirs = {"fastq": fastq_dir, "galaxy": galaxy_dir, "align": align_dir, "work": work_dir,
            "flowcell": fc_dir}
    run_items = add_multiplex_across_lanes(run_info["details"], dirs["fastq"], fc_name)

    # process each flowcell lane
    lane_items = []
    with utils.cpmap(config["algorithm"]["num_cores"]) as cpmap:
        for cur_items in cpmap(process_lane, ((info, fc_name, fc_date, dirs, config)
                                              for info in run_items)):
            lane_items.extend(cur_items)
    with utils.cpmap(config["algorithm"]["num_cores"]) as cpmap:
        for _ in cpmap(process_alignment, lane_items):
            pass
    # process samples, potentially multiplexed across multiple lanes
    sample_files, sample_fastq, sample_info = organize_samples(dirs, fc_name, fc_date, run_items)
    with utils.cpmap(config["algorithm"]["num_cores"]) as cpmap:
        for _ in cpmap(process_sample, ((name, sample_fastq[name], sample_info[name],
                                         bam_files, dirs, config, config_file)
                                        for name, bam_files in sample_files)):
            pass
    write_metrics(run_info, fc_name, fc_date, dirs)

def _get_run_info(fc_name, fc_date, config, run_info_yaml):
    """Retrieve run information from a passed YAML file or the Galaxy API.
    """
    if run_info_yaml and os.path.exists(run_info_yaml):
        log.info("Found YAML samplesheet, using %s instead of Galaxy API" % run_info_yaml)
        with open(run_info_yaml) as in_handle:
            run_details = yaml.load(in_handle)
        return dict(details=run_details, run_id="")
    else:
        log.info("Fetching run details from Galaxy instance")
        galaxy_api = GalaxyApiAccess(config['galaxy_url'], config['galaxy_api_key'])
        return galaxy_api.run_details(fc_name, fc_date)

@utils.map_wrap
def process_lane(info, fc_name, fc_date, dirs, config):
    """Prepare lanes, potentially splitting based on barcodes.
    """
    config = _update_config_w_custom(config, info)

    sample_name = info.get("description", "")
    if (config["algorithm"].get("include_short_name", True) and
            info.get("name", "")):
        sample_name = "%s---%s" % (info.get("name", ""), sample_name)
    genome_build = info.get("genome_build", None)
    multiplex = info.get("multiplex", None)

    log.info("Processing sample: %s; lane %s; reference genome %s; " \
             "researcher %s; analysis method %s" %
             (sample_name, info["lane"], genome_build,
              info.get("researcher", ""), info.get("analysis", "")))
    if multiplex:
        log.debug("Sample %s is multiplexed as: %s" % (sample_name, multiplex))

    full_fastq1, full_fastq2 = get_fastq_files(dirs["fastq"], info['lane'], fc_name)
    lane_name = "%s_%s_%s" % (info['lane'], fc_date, fc_name)
    lane_items = []
    for mname, msample, fastq1, fastq2 in split_by_barcode(full_fastq1,
            full_fastq2, multiplex, lane_name, config):
        mlane_name = "%s_%s" % (lane_name, mname) if mname else lane_name
        if msample is None:
            msample = "%s---%s" % (sample_name, mname)
        lane_items.append((fastq1, fastq2, genome_build, mlane_name, msample,
                           dirs, config))
    return lane_items

@utils.map_wrap
def process_alignment(fastq1, fastq2, genome_build, lane_name, sample, dirs, config):
    """Do an alignment of fastq files, preparing a sorted BAM output file.
    """
    aligner = config["algorithm"].get("aligner", None)
    if os.path.exists(fastq1) and aligner:
        log.info("Aligning lane %s with %s aligner" % (lane_name, aligner))
        align_to_sort_bam(fastq1, fastq2, genome_build, aligner,
                          lane_name, sample, dirs, config)

@utils.map_wrap
def process_sample(sample_name, fastq_files, info, bam_files, dirs, config, config_file):
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
    if config["algorithm"]["recalibrate"]:
        log.info("Recalibrating %s with GATK" % str(sample_name))
        gatk_bam = recalibrate_quality(sort_bam, fastq1, fastq2, sam_ref, config)
        if config["algorithm"]["snpcall"]:
            log.info("SNP genotyping %s with GATK" % str(sample_name))
            vrn_file = run_genotyper(gatk_bam, sam_ref, config)
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

# ## Utility functions

def _get_full_paths(fastq_dir, config, config_file):
    """Retrieve full paths for directories in the case of relative locations.
    """
    fastq_dir = utils.add_full_path(fastq_dir)
    config_dir = utils.add_full_path(os.path.dirname(config_file))
    galaxy_config_file = utils.add_full_path(config["galaxy_config"], config_dir)
    return fastq_dir, os.path.dirname(galaxy_config_file)

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
