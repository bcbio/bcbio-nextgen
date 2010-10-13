#!/usr/bin/env python
"""Perform an automated analysis on a sequencing run using Galaxy information.

Given a directory of solexa output, this retrieves details about the sequencing
run from the Galaxy description, and uses this to perform an initial alignment
and analysis.

Usage:
    automated_initial_analysis.py <YAML config file> <flow cell dir>

Workflow:
    - Retrieve details on a run.
    - Generate fastq files.
    - Align fastq files to reference genome.
    - Generate summary report.

Other items to potentially add:
    - Remap quality scores with GATK.
    - Generate plots of remapped details.

The required elements in the YAML config file are:

galaxy_url: http://galaga/galaxy
galaxy_api_key: your_api_key
galaxy_config: /opt/source/galaxy/web/universe_wsgi.ini
program:
  bowtie: bowtie
  samtools: samtools
  bwa: bwa
  picard: /source/Picard
algorithm:
  aligner: bowtie
  max_errors: 2
  num_cores: 8
  recalibrate: false
"""
import os
import sys
import json
import contextlib
import subprocess
import glob
import copy
import csv
from optparse import OptionParser
from multiprocessing import Pool
import xml.etree.ElementTree as ET

import yaml

from bcbio.solexa.flowcell import (get_flowcell_info, get_fastq_dir)
from bcbio.galaxy.api import GalaxyApiAccess
from bcbio.picard.metrics import PicardMetricsParser
from bcbio.picard import utils

def main(config_file, fc_dir):
    work_dir = os.getcwd()
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    galaxy_api = GalaxyApiAccess(config['galaxy_url'], config['galaxy_api_key'])
    fc_name, fc_date = get_flowcell_info(fc_dir)
    run_info = galaxy_api.run_details(fc_name)
    fastq_dir = get_fastq_dir(fc_dir)
    #print "Generating fastq files"
    #all_lanes = [i['lane'] for i in run_info["details"]]
    #short_fc_name = "%s_%s" % (fc_date, fc_name)
    #fastq_dir = generate_fastq(fc_dir, short_fc_name, all_lanes)
    if config["algorithm"]["num_cores"] > 1:
        pool = Pool(config["algorithm"]["num_cores"])
        try:
            pool.map(_process_wrapper,
                    ((i, fastq_dir, fc_name, fc_date, config, config_file)
                        for i in run_info["details"]))
        except:
            pool.terminate()
            raise
    else:
        map(_process_wrapper,
            ((i, fastq_dir, fc_name, fc_date, config, config_file)
                for i in run_info["details"]))
    write_metrics(run_info, work_dir, fc_dir, fastq_dir)

def process_lane(info, fastq_dir, fc_name, fc_date, config, config_file):
    config = _update_config_w_custom(config, info)
    sample_name = "%s: %s" % (info.get("name", ""),
                              info.get("description", ""))
    genome_build = "%s%s" % (info["genome_build"],
                             config["algorithm"].get("ref_ext", ""))
    print "Processing", info["lane"], genome_build, \
            sample_name, info.get("researcher", ""), \
            info.get("analysis", "")
    lane_name = "%s_%s_%s" % (info['lane'], fc_date, fc_name)
    aligner_to_use = config["algorithm"]["aligner"]
    align_ref, sam_ref = get_genome_ref(genome_build,
            aligner_to_use, os.path.dirname(config["galaxy_config"]))
    fastq1, fastq2 = get_fastq_files(fastq_dir, info['lane'], fc_name)
    print info['lane'], "Aligning with", config["algorithm"]["aligner"]
    if aligner_to_use == "bowtie":
        sam_file = bowtie_to_sam(fastq1, fastq2, align_ref, lane_name,
                config)
    elif aligner_to_use == "bwa":
        sam_file = bwa_align_to_sam(fastq1, fastq2, align_ref, lane_name,
                config)
    elif aligner_to_use == "maq":
        sam_file = maq_align_to_sam(fastq1, fastq2, align_ref, lane_name,
                sample_name, config_file)
    else:
        raise ValueError("Do not recognize aligner: %s" % aligner_to_use)
    print info['lane'], "Converting to sorted BAM file"
    base_bam, sort_bam = sam_to_sort_bam(sam_file, sam_ref, fastq1, fastq2,
            sample_name, config_file)
    bam_to_wig(sort_bam, config, config_file)
    if config["algorithm"]["recalibrate"]:
        print info['lane'], "Recalibrating with GATK"
        dbsnp_file = get_dbsnp_file(config, sam_ref)
        gatk_bam = recalibrate_quality(base_bam, sam_ref,
                dbsnp_file, config["program"]["picard"])
        print info['lane'], "Analyzing recalibration"
        analyze_recalibration(gatk_bam, fastq1, fastq2)
        if config["algorithm"]["snpcall"]:
            print info['lane'], "Providing SNP genotyping with GATK"
            run_genotyper(gatk_bam, sam_ref, dbsnp_file, config_file)

    print info['lane'], "Generating summary files"
    generate_align_summary(sort_bam, fastq1, fastq2, sam_ref,
            config, sample_name, config_file)
    # Cleanup ToDo: 
    # gzip fastq file for storage
    # Remove SAM files

def _process_wrapper(args):
    try:
        return process_lane(*args)
    except KeyboardInterrupt:
        raise Exception

#def generate_fastq(fc_dir, fc_name, all_lanes):
#    fastq_dir = get_fastq_dir(fc_dir)
#    if not fastq_dir == fc_dir and not os.path.exists(fastq_dir):
#        with utils.chdir(os.path.split(fastq_dir)[0]):
#            cl = ["solexa_qseq_to_fastq.py", fc_name,
#                    ",".join(str(l) for l in all_lanes)]
#            subprocess.check_call(cl)
#    return fastq_dir

def bowtie_to_sam(fastq_file, pair_file, ref_file, out_base, config):
    """Before a standard or paired end alignment with bowtie.
    """
    out_file = "%s.sam" % out_base
    if not os.path.exists(out_file):
        cl = [config["program"]["bowtie"]]
        cl += ["-q", "--solexa1.3-quals",
               "-v", config["algorithm"]["max_errors"],
               "-k", 1,
               "-X", 1000, # matches bwa sampe default size
               "-M", 1,
               "--best",
               "--strata",
               "--sam",
               ref_file]
        if pair_file:
            cl += ["-1", fastq_file, "-2", pair_file]
        else:
            cl += [fastq_file]
        cl += [out_file]
        cl = [str(i) for i in cl]
        #print " ".join(cl)
        child = subprocess.Popen(cl)
        child.wait()
    return out_file

def bwa_align_to_sam(fastq_file, pair_file, ref_file, out_base, config):
    """Perform a BWA alignment, generating a SAM file.
    """
    sai1_file = "%s_1.sai" % out_base
    sai2_file = None
    sam_file = "%s.sam" % out_base
    if not os.path.exists(sam_file):
        if not os.path.exists(sai1_file):
            _run_bwa_align(fastq_file, ref_file, sai1_file, config)
        if pair_file:
            sai2_file = "%s_2.sai" % out_base
            if not os.path.exists(sai2_file):
                _run_bwa_align(pair_file, ref_file, sai2_file, config)
        align_type = "sampe" if sai2_file else "samse"
        sam_cl = [config["program"]["bwa"], align_type, ref_file, sai1_file]
        if sai2_file:
            sam_cl.append(sai2_file)
        sam_cl.append(fastq_file)
        if sai2_file:
            sam_cl.append(pair_file)
        with open(sam_file, "w") as out_handle:
            child = subprocess.Popen(sam_cl, stdout=out_handle)
            child.wait()
    return sam_file

def maq_align_to_sam(fastq_file, pair_file, ref_file, out_base,
        sample_name, config_file):
    """Produce a BAM output using Picard to do a maq alignment.
    """
    bam_file = "%s-maq-cal.bam" % out_base
    cl = ["picard_maq_recalibrate.py", "--name=%s" % sample_name,
            config_file, out_base, ref_file, fastq_file]
    if pair_file:
        cl.append(pair_file)
    child = subprocess.Popen(cl)
    child.wait()
    return bam_file

def _run_bwa_align(fastq_file, ref_file, out_file, config):
    aln_cl = [config["program"]["bwa"], "aln",
              "-n %s" % config["algorithm"]["max_errors"],
              "-k %s" % config["algorithm"]["max_errors"],
              ref_file, fastq_file]
    with open(out_file, "w") as out_handle:
        child = subprocess.Popen(aln_cl, stdout=out_handle)
        child.wait()

def sam_to_sort_bam(sam_file, ref_file, fastq1, fastq2, sample_name,
        config_file):
    """Convert SAM file to merged and sorted BAM file.
    """
    base, ext = os.path.splitext(sam_file)
    bam_file = "%s.bam" % base
    sort_bam_file = "%s-sort.bam" % base
    cl = ["picard_sam_to_bam.py", "--name=%s" % sample_name,
            config_file, sam_file, ref_file, fastq1]
    if fastq2:
        cl.append(fastq2)
    subprocess.check_call(cl)
    return bam_file, sort_bam_file

def bam_to_wig(bam_file, config, config_file):
    """Provide a BigWig coverage file of the sorted alignments.
    """
    wig_file = "%s.bigwig" % os.path.splitext(bam_file)[0]
    if not os.path.exists(wig_file):
        cl = [config["analysis"]["towig_script"], bam_file, config_file]
        subprocess.check_call(cl)
    return wig_file

def generate_align_summary(bam_file, fastq1, fastq2, sam_ref, config,
        sample_name, config_file, do_sort=False):
    """Run alignment summarizing script to produce a pdf with align details.
    """
    cl = ["align_summary_report.py", "--name=%s" % sample_name,
            config["program"]["picard"], bam_file, sam_ref, fastq1]
    if fastq2:
        cl.append(fastq2)
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

def recalibrate_quality(bam_file, sam_ref, dbsnp_file, picard_dir):
    """Recalibrate alignments with GATK and provide pdf summary.
    """
    cl = ["picard_gatk_recalibrate.py", picard_dir, sam_ref, bam_file]
    if dbsnp_file:
        cl.append(dbsnp_file)
    subprocess.check_call(cl)
    out_file = glob.glob("%s*gatkrecal.bam" % os.path.splitext(bam_file)[0])[0]
    return out_file

def run_genotyper(bam_file, ref_file, dbsnp_file, config_file):
    """Perform SNP genotyping and analysis using GATK.
    """
    cl = ["gatk_genotyper.py", config_file, ref_file, bam_file]
    if dbsnp_file:
        cl.append(dbsnp_file)
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

def get_fastq_files(directory, lane, fc_name):
    """Retrieve fastq files for the given lane, ready to process.
    """
    files = glob.glob(os.path.join(directory, "%s_*%s*txt*" % (lane, fc_name)))
    files.sort()
    if len(files) > 2 or len(files) == 0:
        raise ValueError("Did not find correct files for %s %s %s %s" %
                (directory, lane, fc_name, files))
    ready_files = []
    for fname in files:
        if fname.endswith(".gz"):
            cl = ["gunzip", fname]
            subprocess.check_call(cl)
            ready_files.append(os.path.splitext(fname)[0])
        else:
            ready_files.append(fname)
    return ready_files[0], (ready_files[1] if len(ready_files) > 1 else None)

def _remap_to_maq(ref_file):
    base_dir = os.path.dirname(os.path.dirname(ref_file))
    name = os.path.basename(ref_file)
    for ext in ["fa", "fasta"]:
        test_file = os.path.join(base_dir, "maq", "%s.%s" % (name, ext))
        if os.path.exists(test_file):
            return test_file
    raise ValueError("Did not find maq file %s" % ref_file)

def get_genome_ref(genome_build, aligner, galaxy_base):
    """Retrieve the reference genome file location from galaxy configuration.
    """
    ref_files = dict(
            bowtie = "bowtie_indices.loc",
            bwa = "bwa_index.loc",
            samtools = "sam_fa_indices.loc",
            maq = "bowtie_indices.loc")
    remap_fns = dict(
            maq = _remap_to_maq
            )
    out_info = []
    for ref_get in [aligner, "samtools"]:
        ref_file = os.path.join(galaxy_base, "tool-data", ref_files[ref_get])
        with open(ref_file) as in_handle:
            for line in in_handle:
                if not line.startswith("#"):
                    parts = line.strip().split()
                    if parts[0] == "index":
                        parts = parts[1:]
                    if parts[0] == genome_build:
                        out_info.append(parts[-1])
                        break
        try:
            out_info[-1] = remap_fns[ref_get](out_info[-1])
        except KeyError:
            pass
        except IndexError:
            raise IndexError("Genome %s not found in %s" % (genome_build,
                ref_file))

    if len(out_info) != 2:
        raise ValueError("Did not find genome reference for %s %s" %
                (genome_build, aligner))
    else:
        return tuple(out_info)

# Output high level summary information for a sequencing run in YAML format
# that can be picked up and loaded into Galaxy.

def write_metrics(run_info, analysis_dir, fc_dir, fastq_dir):
    """Write an output YAML file containing high level sequencing metrics.
    """
    metrics, tab_metrics = summary_metrics(run_info, analysis_dir, fc_dir,
            fastq_dir)
    out_file = os.path.join(analysis_dir, "run_summary.yaml")
    with open(out_file, "w") as out_handle:
        yaml.dump(metrics, out_handle, default_flow_style=False)
    tab_out_file = os.path.join(fc_dir, "run_summary.tsv")
    with open(tab_out_file, "w") as out_handle:
        writer = csv.writer(out_handle, dialect="excel-tab")
        for info in tab_metrics:
            writer.writerow(info)
    return out_file

def summary_metrics(run_info, analysis_dir, fc_dir, fastq_dir):
    """Reformat run and analysis statistics into a YAML-ready format.
    """
    tab_out = []
    out_info = []
    for run in run_info["details"]:
        tab_out.append([run["lane"], run.get("researcher", ""), run["name"],
            run["description"]])
        stats = _metrics_from_stats(_lane_stats(run["lane"], analysis_dir))
        stats.update(_bustard_stats(run["lane"], fastq_dir))
        cur_run_info = dict(
                researcher = run.get("researcher_id", ""),
                sample = run["sample_id"],
                lane = run["lane"],
                request = run_info["run_id"],
                metrics = stats,
                )
        out_info.append(cur_run_info)
    return out_info, tab_out

def _metrics_from_stats(stats):
    """Remap Broad metrics names to our local names.
    """
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

def _bustard_stats(lane_num, fastq_dir):
    """Extract statistics about the flow cell from Bustard outputs.
    """
    sum_file = os.path.join(fastq_dir, os.pardir, "BustardSummary.xml")
    #sum_file = os.path.join(fc_dir, "Data", "Intensities", "BaseCalls",
    #        "BustardSummary.xml")
    with open(sum_file) as in_handle:
        results = ET.parse(in_handle).getroot().find("TileResultsByLane")
        for lane in results:
            if lane.find("laneNumber").text == str(lane_num):
                stats = _collect_cluster_stats(lane)
    return stats

def _collect_cluster_stats(lane):
    """Retrieve total counts on cluster statistics.
    """
    stats = {"Clusters" : 0, "Clusters passed": 0}
    for tile in lane.find("Read").findall("Tile"):
        stats["Clusters"] += int(tile.find("clusterCountRaw").text)
        stats["Clusters passed"] += int(tile.find("clusterCountPF").text)
    return stats

def _lane_stats(lane, work_dir):
    """Parse metrics information from files in the working directory.
    """
    parser = PicardMetricsParser()
    metrics_files = glob.glob(os.path.join(work_dir, "%s_*metrics" % lane))
    metrics = parser.extract_metrics(metrics_files)
    return metrics

# Utility functions

def _update_config_w_custom(config, lane_info):
    """Update the configuration for this lane if a custom analysis is specified.
    """
    config = copy.deepcopy(config)
    custom = config["custom_algorithms"].get(lane_info.get("analysis", None),
            None)
    if custom:
        for key, val in custom.iteritems():
            config["algorithm"][key] = val
    return config

if __name__ == "__main__":
    parser = OptionParser()
    (options, args) = parser.parse_args()
    kwargs = dict()
    main(*args, **kwargs)
