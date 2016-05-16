"""Detect variants with the Platypus Haplotype-Based variant caller.

http://www.well.ox.ac.uk/platypus
https://github.com/andyrimmer/Platypus
"""
import os

import toolz as tz

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import shared as pshared
from bcbio.provenance import do
from bcbio.variation import bamprep, bedutils, vcfutils

def run(align_bams, items, ref_file, assoc_files, region, out_file):
    """Run platypus variant calling, germline whole genome or exome.
    """
    assert out_file.endswith(".vcf.gz")
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            for align_bam in align_bams:
                bam.index(align_bam, items[0]["config"])
            cmd = ["platypus", "callVariants", "--regions=%s" % _subset_regions(region, out_file, items),
                   "--bamFiles=%s" % ",".join(align_bams),
                   "--refFile=%s" % dd.get_ref_file(items[0]), "--output=-",
                   "--logFileName", "/dev/null", "--verbosity=1"]
            resources = config_utils.get_resources("platypus", items[0]["config"])
            if resources.get("options"):
                # normalize options so we can set defaults without overwriting user specified
                for opt in resources["options"]:
                    if "=" in opt:
                        key, val = opt.split("=")
                        cmd.extend([key, val])
                    else:
                        cmd.append(opt)
            # Adjust default filter thresholds to achieve similar sensitivity/specificity to other callers
            # Currently not used after doing more cross validation as they increase false positives
            # which seems to be a major advantage for Platypus users.
            # tuned_opts = ["--hapScoreThreshold", "10", "--scThreshold", "0.99", "--filteredReadsFrac", "0.9",
            #               "--rmsmqThreshold", "20", "--qdThreshold", "0", "--abThreshold", "0.0001",
            #               "--minVarFreq", "0.0", "--assemble", "1"]
            # for okey, oval in utils.partition_all(2, tuned_opts):
            #     if okey not in cmd:
            #         cmd.extend([okey, oval])

            # Avoid filtering duplicates on high depth targeted regions where we don't mark duplicates
            if any(not tz.get_in(["config", "algorithm", "mark_duplicates"], data, True)
                   for data in items):
                cmd += ["--filterDuplicates=0"]
            post_process_cmd = (" | %s | vcfallelicprimitives -t DECOMPOSED --keep-geno | vcffixup - | "
                                "vcfstreamsort | bgzip -c > %s" % (vcfutils.fix_ambiguous_cl(), tx_out_file))
            do.run(" ".join(cmd) + post_process_cmd, "platypus variant calling")
        out_file = vcfutils.bgzip_and_index(out_file, items[0]["config"])
    return out_file

def _subset_regions(region, base_file, items):
    """Subset to a BED file (or genomic region) for calling.
    """
    variant_regions = bedutils.merge_overlaps(bedutils.population_variant_regions(items), items[0])
    target = pshared.subset_variant_regions(variant_regions, region, base_file, items)
    if isinstance(target, basestring) and os.path.isfile(target):
        return target
    else:
        return bamprep.region_to_gatk(target)
