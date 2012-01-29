"""Use Hydra to detect structural variation using discordant read pairs.

Hydra: http://code.google.com/p/hydra-sv/

Pipeline: http://code.google.com/p/hydra-sv/wiki/TypicalWorkflow
"""
import os
import copy
import collections
import subprocess
from contextlib import nested, closing

import pysam
import numpy
from Bio.Seq import Seq

from bcbio import utils, broad
from bcbio.pipeline.alignment import align_to_sort_bam
from bcbio.distributed.transaction import file_transaction

## Prepare alignments to identify discordant pair mappings

def select_unaligned_read_pairs(in_bam, extra, out_dir, config):
    """Retrieve unaligned read pairs from input alignment BAM, as two fastq files.
    """
    runner = broad.runner_from_config(config)
    base, ext = os.path.splitext(os.path.basename(in_bam))
    nomap_bam = os.path.join(out_dir, "{}-{}{}".format(base, extra, ext))
    sort_bam = runner.run_fn("picard_sort", in_bam, "queryname")
    if not utils.file_exists(nomap_bam):
        with file_transaction(nomap_bam) as tx_out:
            runner.run("FilterSamReads", [("INPUT", sort_bam),
                                          ("OUTPUT", tx_out),
                                          ("EXCLUDE_ALIGNED", "true"),
                                          ("WRITE_READS_FILES", "false"),
                                          ("SORT_ORDER", "queryname")])
    has_reads = False
    with closing(pysam.Samfile(nomap_bam, "rb")) as in_pysam:
        for read in in_pysam:
            if read.is_paired:
                has_reads = True
                break
    if has_reads:
        out_fq1, out_fq2 = ["{}-{}.fq".format(os.path.splitext(nomap_bam)[0], i) for i in [1, 2]]
        runner.run_fn("picard_bam_to_fastq", nomap_bam, out_fq1, out_fq2)
        return out_fq1, out_fq2
    else:
        return None, None

def remove_nopairs(in_bam, out_dir):
    """Remove any reads without both pairs present in the file.
    """
    out_bam = os.path.join(out_dir, apply("{}-safepair{}".format,
                                          os.path.splitext(os.path.basename(in_bam))))
    if not utils.file_exists(out_bam):
        read_counts = collections.defaultdict(int)
        with closing(pysam.Samfile(in_bam, "rb")) as in_pysam:
            for read in in_pysam:
                if read.is_paired:
                    read_counts[read.qname] += 1
        with closing(pysam.Samfile(in_bam, "rb")) as in_pysam:
            with file_transaction(out_bam) as tx_out_bam:
                with closing(pysam.Samfile(tx_out_bam, "wb", template=in_pysam)) as out_pysam:
                    for read in in_pysam:
                        if read_counts[read.qname] == 2:
                            out_pysam.write(read)
    return out_bam

def calc_paired_insert_stats(in_bam):
    """Retrieve statistics for paired end read insert distances.
    """
    dists = []
    with closing(pysam.Samfile(in_bam, "rb")) as in_pysam:
        for read in in_pysam:
            if read.is_proper_pair and read.is_read1:
                dists.append(abs(read.isize))
    # remove outliers
    med = numpy.median(dists)
    filter_dists = filter(lambda x: x < med + 10 * med, dists)
    return {"mean": numpy.mean(filter_dists), "std": numpy.std(filter_dists),
            "median": numpy.median(filter_dists)}

def tiered_alignment(in_bam, tier_num, multi_mappers, extra_args,
                     genome_build, pair_stats,
                     work_dir, dirs, config):
    """Perform the alignment of non-mapped reads from previous tier.
    """
    nomap_fq1, nomap_fq2 = select_unaligned_read_pairs(in_bam, "tier{}".format(tier_num),
                                                       work_dir, config)
    if nomap_fq1 is not None:
        base_name = "{}-tier{}out".format(os.path.splitext(os.path.basename(in_bam))[0],
                                          tier_num)
        config = copy.deepcopy(config)
        dirs = copy.deepcopy(dirs)
        config["algorithm"]["multiple_mappers"] = multi_mappers
        config["algorithm"]["extra_align_args"] = ["-i", int(pair_stats["mean"]),
                                               int(pair_stats["std"])] + extra_args
        dirs["align"] = os.path.split(nomap_fq1)[0]
        return align_to_sort_bam(nomap_fq1, nomap_fq2, genome_build, "novoalign",
                                 base_name, base_name,
                                 dirs, config)
    else:
        return None

## Run hydra to identify structural variation breakpoints

@utils.memoize_outfile(".bed")
def convert_bam_to_bed(in_bam, out_file):
    """Convert BAM to bed file using BEDTools.
    """
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            subprocess.check_call(["bamToBed", "-i", in_bam, "-tag", "NM"],
                                  stdout=out_handle)
    return out_file

@utils.memoize_outfile("-pair.bed")
def pair_discordants(in_bed, out_file):
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            subprocess.check_call(["pairDiscordants.py", "-i", in_bed,
                                   "-m", "hydra", "-z", "800"],
                                  stdout=out_handle)
    return out_file

@utils.memoize_outfile("-dedup.bed")
def dedup_discordants(in_bed, out_file):
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            subprocess.check_call(["dedupDiscordants.py", "-i", in_bed, "-s", "3"],
                                  stdout=out_handle)
    return out_file

@utils.memoize_outfile("-hydra.breaks")
def run_hydra(in_bed, pair_stats, out_file):
    with file_transaction(out_file) as tx_out_file:
        subprocess.check_call(["hydra", "-i", in_bed, "-out", tx_out_file,
                               "-mld", str(int(pair_stats["median"])),
                               "-mno", str(int(pair_stats["median"]) +
                                           20 * int(pair_stats["std"]))])
    return out_file

def hydra_breakpoints(in_bam, pair_stats):
    """Detect structural variation breakpoints with hydra.
    """
    in_bed = convert_bam_to_bed(in_bam)
    if os.path.getsize(in_bed) > 0:
        pair_bed = pair_discordants(in_bed)
        dedup_bed = dedup_discordants(pair_bed, pair_stats)
        return run_hydra(dedup_bed, pair_stats)
    else:
        return None

## Top level organizational code

def detect_sv(align_bam, genome_build, dirs, config):
    """Detect structural variation from discordant aligned pairs.
    """
    work_dir = utils.safe_makedir(os.path.join(dirs["work"], "structural"))
    pair_stats = calc_paired_insert_stats(align_bam)
    fix_bam = remove_nopairs(align_bam, work_dir)
    tier2_align = tiered_alignment(fix_bam, "2", True, [],
                                   genome_build, pair_stats,
                                   work_dir, dirs, config)
    if tier2_align:
        tier3_align = tiered_alignment(tier2_align, "3", "Ex 1100", ["-t", "300"],
                                       genome_build, pair_stats,
                                       work_dir, dirs, config)
        if tier3_align:
            hydra_bps = hydra_breakpoints(tier3_align, pair_stats)
