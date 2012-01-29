"""Use Hydra to detect structural variation using discordant read pairs.

Hydra: http://code.google.com/p/hydra-sv/

Pipeline: http://code.google.com/p/hydra-sv/wiki/TypicalWorkflow
"""
import os
import copy
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
    sort_bam = os.path.join(out_dir, "{}-{}{}".format(base, extra, ext))
    runner.run_fn("picard_sort", in_bam, sort_order="queryname", out_file=sort_bam)
    out_fq1, out_fq2 = ["{}-{}.fq".format(os.path.splitext(sort_bam)[0], i) for i in [1, 2]]
    if not utils.file_exists(out_fq1):
        with nested(closing(pysam.Samfile(in_bam, "rb")),
                    open(out_fq1, "wb"), open(out_fq2, "wb")) as (in_pysam, out1, out2):
            for read in in_pysam:
                if read.is_paired and not read.is_proper_pair:
                    cur_out = out1 if read.is_read1 else out2
                    if read.is_reverse:
                        seq = str(Seq(read.seq).reverse_complement())
                        qual = "".join(reversed(read.qual))
                    else:
                        seq = read.seq
                        qual = read.qual
                    cur_out.write("@{name}\n{seq}\n+\n{qual}\n".format(
                        name=read.qname, seq=seq, qual=qual))
    return out_fq1, out_fq2

def calc_paired_insert_stats(in_bam):
    """Retrieve statistics for paired end read insert distances.
    """
    dists = []
    with closing(pysam.Samfile(in_bam, "rb")) as in_pysam:
        for read in in_pysam:
            if read.is_proper_pair and read.is_read1:
                dists.append(read.isize)
    return {"mean": numpy.mean(dists), "std": numpy.std(dists),
            "median": numpy.median(dists)}

def tiered_alignment(in_bam, tier_num, multi_mappers, extra_args,
                     genome_build, pair_stats,
                     work_dir, dirs, config):
    """Perform the alignment of non-mapped reads from previous tier.
    """
    nomap_fq1, nomap_fq2 = select_unaligned_read_pairs(in_bam, "tier{}".format(tier_num),
                                                       work_dir, config)
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
    tier2_align = tiered_alignment(align_bam, "2", True, [],
                                   genome_build, pair_stats,
                                   work_dir, dirs, config)
    tier3_align = tiered_alignment(tier2_align, "3", "Ex 1100", ["-t", "300"],
                                   genome_build, pair_stats,
                                   work_dir, dirs, config)
    hydra_bps = hydra_breakpoints(tier3_align, pair_stats)
