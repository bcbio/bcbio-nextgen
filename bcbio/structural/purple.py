"""PURPLE: Purity and ploidy estimates for somatic tumor/normal samples

https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator
"""
import csv
import os

from bcbio import utils
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
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
    print(het_file)
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

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(dd.get_work_dir(data), "structural",
                                           dd.get_sample_name(data), "purple"))
