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
                           variant_regions=None, known_vrns=None):
    """Generate a list of interval regions for realignment around indels.
    """
    if not known_vrns:
        known_vrns = {}
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
            if known_vrns.get("train_indels"):
                params += ["--known", known_vrns["train_indels"]]
            if deep_coverage:
                params += ["--mismatchFraction", "0.30",
                           "--maxIntervalSize", "650"]
            runner.run_gatk(params, memscale={"direction": "decrease", "magnitude": 2})
    return out_file

def gatk_indel_realignment_cl(runner, align_bam, ref_file, intervals,
                              tmp_dir, region=None, deep_coverage=False,
                              known_vrns=None):
    """Prepare input arguments for GATK indel realignment.
    """
    if not known_vrns:
        known_vrns = {}
    params = ["-T", "IndelRealigner",
              "-I", align_bam,
              "-R", ref_file,
              "-targetIntervals", intervals,
              ]
    if region:
        params += ["-L", region]
    if known_vrns.get("train_indels"):
        params += ["--knownAlleles", known_vrns["train_indels"]]
    if deep_coverage:
        params += ["--maxReadsInMemory", "300000",
                   "--maxReadsForRealignment", str(int(5e5)),
                   "--maxReadsForConsensuses", "500",
                   "--maxConsensuses", "100"]
    return runner.cl_gatk(params, tmp_dir)

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
