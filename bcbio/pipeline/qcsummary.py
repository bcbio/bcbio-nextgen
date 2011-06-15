"""Quality control and summary metrics for next-gen alignments and analysis.
"""
import os
import csv
import copy
import glob
import subprocess
import xml.etree.ElementTree as ET

import yaml

from bcbio.broad.metrics import PicardMetricsParser

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

# Output high level summary information for a sequencing run in YAML format
# that can be picked up and loaded into Galaxy.

def write_metrics(run_info, fc_name, fc_date, dirs):
    """Write an output YAML file containing high level sequencing metrics.
    """
    lane_stats, sample_stats, tab_metrics = summary_metrics(run_info,
            dirs["work"], fc_name, fc_date, dirs["fastq"])
    out_file = os.path.join(dirs["work"], "run_summary.yaml")
    with open(out_file, "w") as out_handle:
        metrics = dict(lanes=lane_stats, samples=sample_stats)
        yaml.dump(metrics, out_handle, default_flow_style=False)
    tab_out_file = os.path.join(dirs["flowcell"], "run_summary.tsv")
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
