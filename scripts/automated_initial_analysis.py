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
import contextlib
import subprocess
import glob
import csv
import copy
import shutil
import collections
from optparse import OptionParser
import xml.etree.ElementTree as ET
import StringIO

import yaml
import logbook

from bcbio.solexa.flowcell import (get_flowcell_info, get_fastq_dir)
from bcbio.galaxy.api import GalaxyApiAccess
from bcbio.broad.metrics import PicardMetricsParser
from bcbio import utils
from bcbio.broad import BroadRunner
from bcbio.log import create_log_handler
from bcbio.pipeline.fastq import get_fastq_files
from bcbio.pipeline.alignment import align_to_sort_bam, get_genome_ref
from bcbio.pipeline.demultiplex import split_by_barcode, add_multiplex_across_lanes

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
    fastq1, fastq2 = _combine_fastq_files(fastq_files, work_dir)
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

def _combine_fastq_files(in_files, work_dir):
    if len(in_files) == 1:
        return in_files[0]
    else:
        cur1, cur2 = in_files[0]
        out1 = os.path.join(work_dir, os.path.basename(cur1))
        out2 = os.path.join(work_dir, os.path.basename(cur2)) if cur2 else None
        if not os.path.exists(out1):
            with open(out1, "a") as out_handle:
                for (cur1, _) in in_files:
                    with open(cur1) as in_handle:
                        shutil.copyfileobj(in_handle, out_handle)
        if out2 and not os.path.exists(out2):
            with open(out2, "a") as out_handle:
                for (_, cur2) in in_files:
                    with open(cur2) as in_handle:
                        shutil.copyfileobj(in_handle, out_handle)
        return out1, out2

def organize_samples(align_dir, fastq_dir, work_dir, fc_name, fc_date, run_items):
    """Organize BAM output files by sample name, handling multiplexing.
    """
    bams_by_sample = collections.defaultdict(list)
    sample_info = dict()
    fastq_by_sample = collections.defaultdict(list)
    for lane_info in run_items:
        multiplex = lane_info.get("multiplex", None)
        if multiplex:
            mfastq_dir = os.path.join(work_dir, "%s_%s_%s_barcode" %
                    (lane_info["lane"], fc_date, fc_name))
            for multi in multiplex:
                name = (lane_info.get("name", ""), lane_info["description"],
                        multi["name"])
                fname = os.path.join(align_dir, "%s_%s_%s_%s-sort.bam" %
                    (lane_info["lane"], fc_date, fc_name, multi["barcode_id"]))
                if os.path.exists(fname):
                    bams_by_sample[name].append(fname)
                    sample_info[name] = lane_info
                    fastq_by_sample[name].append(get_fastq_files(mfastq_dir,
                        lane_info["lane"], fc_name, multi["barcode_id"]))
        else:
            name = (lane_info.get("name", ""), lane_info["description"])
            fname = os.path.join(align_dir, "%s_%s_%s-sort.bam" %
                    (lane_info["lane"], fc_date, fc_name))
            if os.path.exists(fname):
                bams_by_sample[name].append(fname)
                sample_info[name] = lane_info
                fastq_by_sample[name].append(get_fastq_files(fastq_dir,
                    lane_info["lane"], fc_name))
    return sorted(bams_by_sample.items()), dict(fastq_by_sample), sample_info

def merge_bam_files(bam_files, work_dir, config):
    """Merge multiple BAM files from a sample into a single BAM for processing.
    """
    out_file = os.path.join(work_dir, os.path.basename(bam_files[0]))
    picard = BroadRunner(config["program"]["picard"],
                         max_memory=config["algorithm"].get("java_memory", ""))
    picard.run_fn("picard_merge", bam_files, out_file)
    for b in bam_files:
        utils.save_diskspace(b, "BAM merged to %s" % out_file, config)
    return out_file

def bam_to_wig(bam_file, config, config_file):
    """Provide a BigWig coverage file of the sorted alignments.
    """
    wig_file = "%s.bigwig" % os.path.splitext(bam_file)[0]
    if not (os.path.exists(wig_file) and os.path.getsize(wig_file) > 0):
        cl = [config["analysis"]["towig_script"], bam_file, config_file]
        subprocess.check_call(cl)
    return wig_file

def generate_align_summary(bam_file, is_paired, sam_ref, config,
        sample_name, config_file, do_sort=False):
    """Run alignment summarizing script to produce a pdf with align details.
    """
    sample_name = " : ".join(sample_name)
    cl = ["align_summary_report.py", "--name=%s" % sample_name,
            config["program"]["picard"], bam_file, sam_ref]
    if is_paired:
        cl.append("--paired")
    bait = config["algorithm"].get("hybrid_bait", "")
    target = config["algorithm"].get("hybrid_target", "")
    if bait and target:
        base_dir = os.path.dirname(os.path.dirname(sam_ref))
        cl.append("--bait=%s" % os.path.join(base_dir, bait))
        cl.append("--target=%s" % os.path.join(base_dir, target))
    if do_sort:
        cl.append("--sort")
    cl.append("--config=%s" % config_file)
    subprocess.check_call(cl)

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

# Output high level summary information for a sequencing run in YAML format
# that can be picked up and loaded into Galaxy.

def write_metrics(run_info, analysis_dir, fc_dir, fc_name, fc_date,
        fastq_dir):
    """Write an output YAML file containing high level sequencing metrics.
    """
    lane_stats, sample_stats, tab_metrics = summary_metrics(run_info,
            analysis_dir, fc_name, fc_date, fastq_dir)
    out_file = os.path.join(analysis_dir, "run_summary.yaml")
    with open(out_file, "w") as out_handle:
        metrics = dict(lanes=lane_stats, samples=sample_stats)
        yaml.dump(metrics, out_handle, default_flow_style=False)
    tab_out_file = os.path.join(fc_dir, "run_summary.tsv")
    try:
        with open(tab_out_file, "w") as out_handle:
            writer = csv.writer(out_handle, dialect="excel-tab")
            for info in tab_metrics:
                writer.writerow(info)
    # If on NFS mounted directory can fail due to filesystem or permissions
    # errors. That's okay, we'll just not write the file.
    except IOError:
        pass
    return out_file

def summary_metrics(run_info, analysis_dir, fc_name, fc_date, fastq_dir):
    """Reformat run and analysis statistics into a YAML-ready format.
    """
    tab_out = []
    lane_info = []
    sample_info = []
    for run in run_info["details"]:
        tab_out.append([run["lane"], run.get("researcher", ""),
            run.get("name", ""), run.get("description")])
        base_info = dict(
                researcher = run.get("researcher_id", ""),
                sample = run.get("sample_id", ""),
                lane = run["lane"],
                request = run_info["run_id"])
        cur_lane_info = copy.deepcopy(base_info)
        cur_lane_info["metrics"] = _bustard_stats(run["lane"], fastq_dir,
                fc_date)
        lane_info.append(cur_lane_info)
        for barcode in run.get("multiplex", [None]):
            cur_name = "%s_%s_%s" % (run["lane"], fc_date, fc_name)
            if barcode:
                cur_name = "%s_%s-" % (cur_name, barcode["barcode_id"])
            stats = _metrics_from_stats(_lane_stats(cur_name, analysis_dir))
            if stats:
                cur_run_info = copy.deepcopy(base_info)
                cur_run_info["metrics"] = stats
                cur_run_info["barcode_id"] = str(barcode["barcode_id"]) if barcode else ""
                cur_run_info["barcode_type"] = (str(barcode.get("barcode_type", ""))
                                                if barcode else "")
                sample_info.append(cur_run_info)
    return lane_info, sample_info, tab_out

def _metrics_from_stats(stats):
    """Remap Broad metrics names to our local names.
    """
    if stats:
        s_to_m = dict(
                AL_MEAN_READ_LENGTH = 'Read length',
                AL_TOTAL_READS = 'Reads',
                AL_PF_READS_ALIGNED = 'Aligned',
                DUP_READ_PAIR_DUPLICATES = 'Pair duplicates'
                )
        metrics = dict()
        for stat_name, metric_name in s_to_m.iteritems():
            metrics[metric_name] = stats[stat_name]
        return metrics

def _bustard_stats(lane_num, fastq_dir, fc_date):
    """Extract statistics about the flow cell from Bustard outputs.
    """
    sum_file = os.path.join(fastq_dir, os.pardir, "BustardSummary.xml")
    #sum_file = os.path.join(fc_dir, "Data", "Intensities", "BaseCalls",
    #        "BustardSummary.xml")
    stats = dict()
    if os.path.exists(sum_file):
        with open(sum_file) as in_handle:
            results = ET.parse(in_handle).getroot().find("TileResultsByLane")
            for lane in results:
                if lane.find("laneNumber").text == str(lane_num):
                    stats = _collect_cluster_stats(lane)
    read_stats = _calc_fastq_stats(fastq_dir, lane_num, fc_date)
    stats.update(read_stats)
    return stats

def _calc_fastq_stats(fastq_dir, lane_num, fc_date):
    """Grab read length from fastq; could provide distribution if non-equal.
    """
    stats = dict()
    fastq_files = glob.glob(os.path.join(fastq_dir, "%s_%s*" % (lane_num,
        fc_date)))
    if len(fastq_files) > 0:
        fastq_file = sorted(fastq_files)[-1]
        with open(fastq_file) as in_handle:
            line = in_handle.readline()
            stats["Read length"] = len(in_handle.readline().strip())
    return stats

def _collect_cluster_stats(lane):
    """Retrieve total counts on cluster statistics.
    """
    stats = {"Clusters" : 0, "Clusters passed": 0}
    for tile in lane.find("Read").findall("Tile"):
        stats["Clusters"] += int(tile.find("clusterCountRaw").text)
        stats["Clusters passed"] += int(tile.find("clusterCountPF").text)
    return stats

def _lane_stats(cur_name, work_dir):
    """Parse metrics information from files in the working directory.
    """
    parser = PicardMetricsParser()
    metrics_files = glob.glob(os.path.join(work_dir, "%s*metrics" % cur_name))
    metrics = parser.extract_metrics(metrics_files)
    return metrics

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
