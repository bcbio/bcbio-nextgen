"""Perform realignment of BAM files around indels using the GATK toolkit.
"""
import os
import shutil
from contextlib import closing

import pysam

from bcbio import bam, broad
from bcbio.bam import ref
from bcbio.log import logger
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.pipeline.shared import subset_bam_by_region, subset_variant_regions
from bcbio.provenance import do

# ## GATK realignment

def gatk_realigner_targets(runner, align_bam, ref_file, config, dbsnp=None,
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
        with file_transaction(config, out_file) as tx_out_file:
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
            runner.run_gatk(params, memscale={"direction": "decrease", "magnitude": 2})
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
                           region=None, out_file=None, deep_coverage=False,
                           config=None):
    """Perform realignment of BAM file in specified regions
    """
    if out_file is None:
        out_file = "%s-realign.bam" % os.path.splitext(align_bam)[0]
    if not file_exists(out_file):
        with tx_tmpdir(config) as tmp_dir:
            with file_transaction(config, out_file) as tx_out_file:
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
    bam.index(align_bam, config)
    runner.run_fn("picard_index_ref", ref_file)
    ref.fasta_idx(ref_file)
    if region:
        align_bam = subset_bam_by_region(align_bam, region, config, out_file)
        bam.index(align_bam, config)
    if has_aligned_reads(align_bam, region):
        variant_regions = config["algorithm"].get("variant_regions", None)
        realign_target_file = gatk_realigner_targets(runner, align_bam,
                                                     ref_file, config, dbsnp, region,
                                                     out_file, deep_coverage,
                                                     variant_regions)
        realign_bam = gatk_indel_realignment(runner, align_bam, ref_file,
                                             realign_target_file, region,
                                             out_file, deep_coverage, config=config)
        # No longer required in recent GATK (> Feb 2011) -- now done on the fly
        # realign_sort_bam = runner.run_fn("picard_fixmate", realign_bam)
        return realign_bam
    elif out_file:
        shutil.copy(align_bam, out_file)
        return out_file
    else:
        return align_bam

# ## Utilities

def has_aligned_reads(align_bam, region=None):
    """Check if the aligned BAM file has any reads in the region.

    region can be a chromosome string ("chr22"),
    a tuple region (("chr22", 1, 100)) or a file of regions.
    """
    import pybedtools
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
                        return True
                else:
                    for item in cur_bam.fetch(region[0], int(region[1]), int(region[2])):
                        return True
        else:
            for item in cur_bam:
                if not item.is_unmapped:
                    return True
    return False
