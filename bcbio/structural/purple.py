"""PURPLE: Purity and ploidy estimates for somatic tumor/normal samples

https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator
"""
import csv
import os

import toolz as tz

from bcbio import utils
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.distributed.transaction import file_transaction
from bcbio.variation import vcfutils

def run(items):
    paired = vcfutils.get_paired(items)
    if not paired or not paired.normal_name:
        logger.info("Skipping PURPLE; need tumor/normal somatic calls in batch: %s" %
                    " ".join([dd.get_sample_name(d) for d in items]))
        return items
    work_dir = _sv_workdir(paired.tumor_data)
    from bcbio import heterogeneity
    het_file = _amber_het_file(heterogeneity.get_variants(paired.tumor_data), work_dir, paired)
    depth_file = _run_cobalt(paired, work_dir)
    print(het_file, depth_file)
    return items

class OutWriter:
    def __init__(self, out_handle):
        self.writer = csv.writer(out_handle, dialect="excel-tab")

    def write_header(self):
        self.writer.writerow(["Chromosome", "Position", "TumorBAF", "TumorModifiedBAF", "TumorDepth",
                                "NormalBAF", "NormalModifiedBAF", "NormalDepth"])

    def _normalize_baf(self, baf):
        """Provide normalized BAF in the same manner as Amber, relative to het.

        https://github.com/hartwigmedical/hmftools/blob/637e3db1a1a995f4daefe2d0a1511a5bdadbeb05/hmf-common/src/main/java/com/hartwig/hmftools/common/amber/AmberBAF.java#L16
        """
        return 0.5 + abs(baf - 0.5)

    def write_row(self, rec, stats):
        self.writer.writerow([rec.chrom, rec.pos,
                              stats["tumor"]["freq"], self._normalize_baf(stats["tumor"]["freq"]),
                              stats["tumor"]["depth"],
                              stats["normal"]["freq"], self._normalize_baf(stats["normal"]["freq"]),
                              stats["normal"]["depth"]])

def _amber_het_file(vrn_files, work_dir, paired):
    """Create file of BAFs in normal heterozygous positions compatible with AMBER.

    https://github.com/hartwigmedical/hmftools/tree/master/amber
    https://github.com/hartwigmedical/hmftools/blob/637e3db1a1a995f4daefe2d0a1511a5bdadbeb05/hmf-common/src/test/resources/amber/new.amber.baf
    """
    assert vrn_files, "Did not find compatible variant calling files for TitanCNA inputs"
    from bcbio.heterogeneity import bubbletree

    prep_file = bubbletree.prep_vrn_file(vrn_files[0]["vrn_file"], vrn_files[0]["variantcaller"],
                                         work_dir, paired, OutWriter)
    amber_dir = utils.safe_makedir(os.path.join(work_dir, "amber"))
    out_file = os.path.join(amber_dir, "%s.amber.baf" % dd.get_sample_name(paired.tumor_data))
    utils.symlink_plus(prep_file, out_file)
    return out_file

def _run_cobalt(paired, work_dir):
    """Run Cobalt for counting read depth across genomic windows.

    PURPLE requires even 1000bp windows so use integrated counting solution
    directly rather than converting from CNVkit calculations. If this approach
    is useful should be moved upstream to be available to other tools as
    an input comparison.

    https://github.com/hartwigmedical/hmftools/tree/master/count-bam-lines
    """
    pass

def _cobalt_ratio_file(paired, work_dir):
    """Convert CNVkit binning counts into cobalt ratio output.

    This contains read counts plus normalization for GC, from section 7.2
    "Determine read depth ratios for tumor and reference genomes"

    https://www.biorxiv.org/content/biorxiv/early/2018/09/20/415133.full.pdf

    Since CNVkit cnr files already have GC bias correction, we re-center
    the existing log2 ratios to be around 1, rather than zero, which matches
    the cobalt expectations.

    XXX This doesn't appear to be a worthwhile direction since PURPLE requires
    1000bp even binning. We'll leave this here as a starting point for future
    work but work on using cobalt directly.
    """
    cobalt_dir = utils.safe_makedir(os.path.join(work_dir, "cobalt"))
    out_file = os.path.join(cobalt_dir, "%s.cobalt" % dd.get_sample_name(paired.tumor_data))
    if not utils.file_exists(out_file):
        cnr_file = tz.get_in(["depth", "bins", "normalized"], paired.tumor_data)
        with file_transaction(paired.tumor_data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                writer = csv.writer(out_handle, delimiter="\t")
                writer.writerow(["Chromosome", "Position", "ReferenceReadCount", "TumorReadCount",
                                 "ReferenceGCRatio", "TumorGCRatio", "ReferenceGCDiploidRatio"])
        print(cnr_file)
        raise NotImplementedError
    return out_file

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(dd.get_work_dir(data), "structural",
                                           dd.get_sample_name(data), "purple"))
