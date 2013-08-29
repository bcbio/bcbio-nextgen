"""Perform realignment of BAM files around indels using the GATK toolkit.
"""
import os
import shutil
from contextlib import closing

import pybedtools
import pysam

from bcbio import broad
from bcbio.log import logger
from bcbio.utils import curdir_tmpdir, file_exists, save_diskspace
from bcbio.distributed.transaction import file_transaction
from bcbio.distributed.split import parallel_split_combine
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import (process_bam_by_chromosome, configured_ref_file,
                                   write_nochr_reads, subset_bam_by_region,
                                   subset_variant_regions)
from bcbio.provenance import do

# ## gkno Marth lab realignment

def gkno_realigner_cl(ref_file, config):
    """Prepare commandline for Marth lab realignment tools.
    Assumes feeding to piped input and output so doesn't not manage
    readying or writing from disk.
    """
    ogap = config_utils.get_program("ogap", config)
    bamleftalign = config_utils.get_program("bamleftalign", config)
    cmd = ("{ogap} --repeat-gap-extend 25 --soft-clip-qsum 20 "
           "  --fasta-reference {ref_file} --entropy-gap-open "
           "  --mismatch-qsum 20 --soft-clip-limit 0 "
           "| {bamleftalign} --fasta-reference {ref_file} ")
    return cmd.format(**locals())

def gkno_realigner(align_bam, ref_file, config, dbsnp=None, region=None,
                   out_file=None, deep_coverage=False):
    """Perform realignment using commandline tools from the Marth lab.

    Runs bamtools filter -> ogap -> bamleftalign

    http://blog.gkno.me/post/32258606906/call-short-variants
    """
    if not out_file:
        base, ext = os.path.splitext(align_bam)
        out_file = "%s-realign%s%s" % (base, ("-%s" % region if region else ""), ext)
    bamtools = config_utils.get_program("bamtools", config)
    realign_cmd = gkno_realigner_cl(ref_file, config)
    region = "-region %s" % region if region else ""

    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            cmd = ("{bamtools} filter -in {align_bam} {region} "
                   "| {realign_cmd} > {tx_out_file}")
            do.run(cmd.format(**locals()), "gkno realignment", {})
    return out_file

# ## GATK realignment

def gatk_realigner_targets(runner, align_bam, ref_file, dbsnp=None,
                           region=None, out_file=None, deep_coverage=False,
                           variant_regions=None):
    """Generate a list of interval regions for realignment around indels.
    """
    if out_file:
        out_file = "%s.intervals" % os.path.splitext(out_file)[0]
    else:
        out_file = "%s-realign.intervals" % os.path.splitext(align_bam)[0]
    # check only for file existence; interval files can be empty after running
    # on small chromosomes, so don't rerun in those cases
    if not os.path.exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            logger.debug("GATK RealignerTargetCreator: %s %s" %
                         (os.path.basename(align_bam), region))
            params = ["-T", "RealignerTargetCreator",
                      "-I", align_bam,
                      "-R", ref_file,
                      "-o", tx_out_file,
                      "-l", "INFO",
                      ]
            region = subset_variant_regions(variant_regions, region, tx_out_file)
            if region:
                params += ["-L", region, "--interval_set_rule", "INTERSECTION"]
            if dbsnp:
                params += ["--known", dbsnp]
            if deep_coverage:
                params += ["--mismatchFraction", "0.30",
                           "--maxIntervalSize", "650"]
            runner.run_gatk(params)
    return out_file

def gatk_indel_realignment_cl(runner, align_bam, ref_file, intervals,
                              tmp_dir, region=None, deep_coverage=False):
    """Prepare input arguments for GATK indel realignment.
    """
    params = ["-T", "IndelRealigner",
              "-I", align_bam,
              "-R", ref_file,
              "-targetIntervals", intervals,
              ]
    if region:
        params += ["-L", region]
    if deep_coverage:
        params += ["--maxReadsInMemory", "300000",
                   "--maxReadsForRealignment", str(int(5e5)),
                   "--maxReadsForConsensuses", "500",
                   "--maxConsensuses", "100"]
    return runner.cl_gatk(params, tmp_dir)

def gatk_indel_realignment(runner, align_bam, ref_file, intervals,
                           region=None, out_file=None, deep_coverage=False):
    """Perform realignment of BAM file in specified regions
    """
    if out_file is None:
        out_file = "%s-realign.bam" % os.path.splitext(align_bam)[0]
    if not file_exists(out_file):
        with curdir_tmpdir() as tmp_dir:
            with file_transaction(out_file) as tx_out_file:
                logger.info("GATK IndelRealigner: %s %s" %
                            (os.path.basename(align_bam), region))
                cl = gatk_indel_realignment_cl(runner, align_bam, ref_file, intervals,
                                                   tmp_dir, region, deep_coverage)
                cl += ["-o", tx_out_file]
                do.run(cl, "GATK indel realignment", {})
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
    if region:
        align_bam = subset_bam_by_region(align_bam, region, out_file)
        runner.run_fn("picard_index", align_bam)
    if has_aligned_reads(align_bam, region):
        variant_regions = config["algorithm"].get("variant_regions", None)
        realign_target_file = gatk_realigner_targets(runner, align_bam,
                                                     ref_file, dbsnp, region,
                                                     out_file, deep_coverage,
                                                     variant_regions)
        realign_bam = gatk_indel_realignment(runner, align_bam, ref_file,
                                             realign_target_file, region,
                                             out_file, deep_coverage)
        # No longer required in recent GATK (> Feb 2011) -- now done on the fly
        # realign_sort_bam = runner.run_fn("picard_fixmate", realign_bam)
        return realign_bam
    elif out_file:
        shutil.copy(align_bam, out_file)
        return out_file
    else:
        return align_bam

# ## High level functionality to run realignments in parallel

def has_aligned_reads(align_bam, region=None):
    """Check if the aligned BAM file has any reads in the region.

    region can be a chromosome string ("chr22"),
    a tuple region (("chr22", 1, 100)) or a file of regions.
    """
    has_items = False
    if region is not None:
        if isinstance(region, basestring) and os.path.isfile(region):
            regions = [tuple(r) for r in pybedtools.BedTool(region)]
        else:
            regions = [region]
    with closing(pysam.Samfile(align_bam, "rb")) as cur_bam:
        if region is not None:
            for region in regions:
                if isinstance(region, basestring):
                    for item in cur_bam.fetch(region):
                        has_items = True
                        break
                else:
                    for item in cur_bam.fetch(region[0], int(region[1]), int(region[2])):
                        has_items = True
                        break
        else:
            for item in cur_bam:
                if not item.is_unmapped:
                    has_items = True
                    break
    return has_items

def parallel_realign_sample(sample_info, parallel_fn):
    """Realign samples, running in parallel over individual chromosomes.
    """
    to_process = []
    finished = []
    for x in sample_info:
        if x[0]["config"]["algorithm"].get("realign", True):
            to_process.append(x)
        else:
            finished.append(x)
    if len(to_process) > 0:
        file_key = "work_bam"
        split_fn = process_bam_by_chromosome("-realign.bam", file_key,
                                           default_targets=["nochr"])
        processed = parallel_split_combine(to_process, split_fn, parallel_fn,
                                           "realign_sample", "combine_bam",
                                           file_key, ["config"])
        finished.extend(processed)
    return finished

_realign_approaches = {"gatk": gatk_realigner,
                       "gkno": gkno_realigner}

def realign_sample(data, region=None, out_file=None):
    """Realign sample BAM file at indels.
    """
    realigner = data["config"]["algorithm"].get("realign", True)
    realigner = "gatk" if realigner is True else realigner
    realign_fn = _realign_approaches[realigner] if realigner else None

    if realign_fn:
        logger.info("Realigning %s with %s: %s %s" % (data["name"], realigner,
                                                      os.path.basename(data["work_bam"]),
                                                      region))
        sam_ref = data["sam_ref"]
        config = data["config"]
        if region == "nochr":
            realign_bam = write_nochr_reads(data["work_bam"], out_file)
        else:
            realign_bam = realign_fn(data["work_bam"], sam_ref, config,
                                     configured_ref_file("dbsnp", config, sam_ref),
                                     region, out_file)
        if region is None:
            save_diskspace(data["work_bam"], "Realigned to %s" % realign_bam,
                           config)
        data["work_bam"] = realign_bam
    return [data]
