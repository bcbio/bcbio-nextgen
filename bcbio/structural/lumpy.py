"""Structural variation detection for split and paired reads using lumpy.

https://github.com/arq5x/lumpy-sv#paired-end-and-split-read-alignment-using-bwa-mem
"""
import contextlib
import itertools
import os
import sys
import random

import pybedtools
import pysam
import yaml

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.structural import hydra

# ## Paired end preparation

def _extract_pe_aligns(in_bam, work_dir, data):
    """Extract paired end alignments from input BAM file.
    """
    out_file = os.path.join(work_dir, "%s-pe%s" % os.path.splitext(os.path.basename(in_bam)))
    samtools = config_utils.get_program("samtools", data["config"])
    #resources = config_utils.get_resources("samtools", data["config"])
    #num_cores = data["config"]["algorithm"].get("num_cores", 1)
    #max_mem = resources.get("memory", "1G")
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            cmd = ("{samtools} view -u -F 0x0100 {in_bam} | "
                   "{samtools} view -u -F 0x0004 - | "
                   "{samtools} view -u -F 0x0008 - | "
                   "{samtools} view -b -F 0x0400 - > {tx_out_file} ")
            do.run(cmd.format(**locals()), "Lumpy prep: paired end", data)
    return out_file

def _random_region(region_bed):
    """Generator of "random" genomic regions from the available list.

    Uses a consistent seed to ensure same results for identical sets of callable regions.
    """
    random.seed(42)
    regions = [tuple(r) for r in pybedtools.BedTool(region_bed)]
    while True:
        chrom, start, end = random.choice(regions)
        position = random.randint(int(start), int(end))
        yield (chrom, position)

def _sample_insert_sizes(n, in_bam, callable_regions):
    """Randomly sample insert sizes from the input BAM file in callable regions.
    """
    out = []
    region_gen = _random_region(callable_regions)
    with contextlib.closing(pysam.Samfile(in_bam, "rb")) as in_pysam:
        while len(out) < n:
            chrom, start = region_gen.next()
            for read in in_pysam.fetch(chrom, start, start+100):
                if read.is_proper_pair and read.is_read1:
                    out.append(abs(read.tlen))
                    break
    return out

def _pe_insert_hist(inserts, stats, X, read_length, in_bam):
    """Create an insert size histogram to pass to lumpy.
    """
    start = read_length
    end = int(stats["mean"] + X * stats["std"])
    vals = [0] * (end - start + 1)
    for insert in inserts:
        if insert >= start and insert <= end:
            vals[insert-start] += 1
    total = float(sum(vals))
    return [x/total for x in vals]

def _pe_stats(in_bam, data, read_length, z):
    """Calculate paired end stats for running lumpy.
    """
    n = 1000 # number of reads to estimate insert size
    out_file = "%s.histo" % os.path.splitext(in_bam)[0]
    stats_file = "%s.stats" % os.path.splitext(in_bam)[0]
    if not utils.file_exists(out_file) or not utils.file_exists(stats_file):
        bam.index(in_bam, data["config"])
        inserts = _sample_insert_sizes(n, in_bam, data["config"]["algorithm"]["callable_regions"])
        stats = hydra.insert_size_stats(inserts)
        with file_transaction(out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                for i, val in enumerate(_pe_insert_hist(inserts, stats, z, read_length, in_bam)):
                    out_handle.write("%s\t%s\n" % (i, val))
        with file_transaction(stats_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                yaml.dump(stats, out_handle, default_flow_style=False, allow_unicode=False)
    with open(stats_file) as in_handle:
        stats = yaml.load(in_handle)
    return out_file, stats

def bam_read_length(in_bam):
    """Calculate max read length from first 500 reads.
    """
    with contextlib.closing(pysam.Samfile(in_bam, "rb")) as in_pysam:
        return max(x.qlen for x in itertools.islice(in_pysam, 500))

# ## Split read preparation

def _extract_sr_aligns(in_bam, work_dir, data):
    """Retrieve split-read alignments from input BAM file.
    """
    out_file = os.path.join(work_dir, "%s-sr%s" % os.path.splitext(os.path.basename(in_bam)))
    samtools = config_utils.get_program("samtools", data["config"])
    lumpy = config_utils.get_program("lumpy", data["config"])
    extract_sr_script = os.path.join(os.path.dirname(lumpy), os.pardir, "share",
                                     "lumpy-sv", "extractSplitReads_BwaMem")
    python = sys.executable
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            cmd = ("{samtools} view -h {in_bam} | "
                   "{python} {extract_sr_script} -i stdin | "
                   "{samtools} view -Sb - > {tx_out_file}")
            do.run(cmd.format(**locals()), "Lumpy prep: split read", data)
    return out_file

# ## Lumpy main

def _run_lumpy(pe_data, sr_data, work_dir, data):
    """Run lumpy-sv
    """
    out_file = os.path.join(work_dir, "%s-pesr.bedpe"
                            % os.path.splitext(os.path.basename(data["work_bam"]))[0])
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            lumpy = config_utils.get_program("lumpy", data["config"])
            pe_data_str = ",".join("%s:%s" % (k, v) for k, v in pe_data.iteritems())
            sr_data_str = ",".join("%s:%s" % (k, v) for k, v in sr_data.iteritems())
            cmd = ("{lumpy} -mw 4 -tt 1e-3 -pe {pe_data_str} -sr {sr_data_str} > {tx_out_file}")
            do.run(cmd.format(**locals()), "Lumpy", data)
    return out_file

def run(data):
    """Perform detection of structural variations with lumpy, using bwa-mem alignment.
    """
    z = 4 # range of std deviations considered normal
    if data["config"]["algorithm"].get("aligner") != "bwa":
        raise ValueError("Require bwa-mem alignment input for lumpy structural variation detection")
    work_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural", data["name"][-1],
                                               "lumpy"))
    pe_bam = _extract_pe_aligns(data["work_bam"], work_dir, data)
    read_length = bam_read_length(pe_bam)
    pe_histo, pe_stats = _pe_stats(pe_bam, data, read_length, z)
    sr_bam = _extract_sr_aligns(data["work_bam"], work_dir, data)
    pe_data = {"bam_file": pe_bam, "histo_file": pe_histo, "mean": pe_stats["mean"],
               "stdev": pe_stats["std"], "read_length": read_length,
               "min_non_overlap": read_length, "discordant_z": z,
               "back_distance": 20, "weight": 1, "id": 1, "min_mapping_threshold": 1}
    sr_data = {"bam_file": sr_bam, "back_distance": 20, "weight": 1, "id": 1}
    pebed_file = _run_lumpy(pe_data, sr_data, work_dir, data)
    return {"bed": pebed_file}
