"""Perform realignment of BAM files around indels using the GATK toolkit.
"""
import os

import pysam

from bcbio import broad
from bcbio.pipeline import log
from bcbio.utils import curdir_tmpdir, file_transaction
from bcbio.distributed.split import parallel_split_combine
from bcbio.pipeline.shared import (split_bam_by_chromosome, ref_genome_info,
                                   configured_ref_file)

# ## Realignment runners with GATK specific arguments

def gatk_realigner_targets(runner, align_bam, ref_file, dbsnp=None,
                           region=None, out_file=None, deep_coverage=False):
    """Generate a list of interval regions for realignment around indels.
    """
    if out_file:
        out_file = "%s.intervals" % os.path.splitext(out_file)[0]
    else:
        out_file = "%s-realign.intervals" % os.path.splitext(align_bam)[0]
    params = ["-T", "RealignerTargetCreator",
              "-I", align_bam,
              "-R", ref_file,
              "-o", out_file,
              "-l", "INFO",
              ]
    if region:
        params += ["-L", region]
    if dbsnp:
        params += ["--known", dbsnp]
    if deep_coverage:
        params += ["--mismatchFraction", "0.30",
                   "--maxIntervalSize", "650"]
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        with file_transaction(out_file):
            runner.run_gatk(params)
    return out_file

def gatk_indel_realignment(runner, align_bam, ref_file, intervals,
                           region=None, out_file=None, deep_coverage=False):
    """Perform realignment of BAM file in specified regions
    """
    if out_file is None:
        out_file = "%s-realign.bam" % os.path.splitext(align_bam)[0]
    params = ["-T", "IndelRealigner",
              "-I", align_bam,
              "-R", ref_file,
              "-targetIntervals", intervals,
              "-o", out_file,
              "-l", "INFO",
              ]
    if region:
        params += ["-L", region]
    if deep_coverage:
        params += ["--maxReadsInMemory", "300000",
                   "--maxReadsForRealignment", str(int(5e5)),
                   "--maxReadsForConsensuses", "500",
                   "--maxConsensuses", "100"]
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        with curdir_tmpdir() as tmp_dir:
            with file_transaction(out_file):
                runner.run_gatk(params, tmp_dir)
    return out_file

def gatk_realigner(align_bam, ref_file, config, dbsnp=None, region=None,
                   out_file=None, deep_coverage=False):
    """Realign a BAM file around indels using GATK, returning sorted BAM.
    """
    runner = broad.runner_from_config(config)
    runner.run_fn("picard_index", align_bam)
    runner.run_fn("picard_index_ref", ref_file)
    if not os.path.exists("%s.fai" % ref_file):
        pysam.faidx(ref_file)
    realign_target_file = gatk_realigner_targets(runner, align_bam,
                                                 ref_file, dbsnp, region,
                                                 out_file, deep_coverage)
    realign_bam = gatk_indel_realignment(runner, align_bam, ref_file,
                                         realign_target_file, region,
                                         out_file, deep_coverage)
    # No longer required in recent GATK (> Feb 2011) -- now done on the fly
    # realign_sort_bam = runner.run_fn("picard_fixmate", realign_bam)
    return realign_bam

# ## High level functionality to run realignments in parallel

def parallel_realign_sample(sample_info, parallel_fn):
    """Realign samples, running in parallel over individual chromosomes.
    """
    data = sample_info[0]
    if data["config"]["algorithm"]["snpcall"]:
        file_key = "work_bam"
        split_fn = split_bam_by_chromosome("-realign.bam", file_key)
        return parallel_split_combine(sample_info, split_fn, parallel_fn,
                                      "realign_sample", "combine_bam",
                                      file_key, ["config"])
    else:
        return sample_info

def realign_sample(data, region=None, out_file=None):
    """Realign sample BAM file at indels.
    """
    log.info("Realigning %s with GATK" % str(data["name"]))
    if data["config"]["algorithm"]["snpcall"]:
        sam_ref = data["sam_ref"]
        config = data["config"]
        data["work_bam"] = gatk_realigner(data["work_bam"], sam_ref, config,
                                          configured_ref_file("dbsnp", config, sam_ref),
                                          region, out_file)
    return [data]
