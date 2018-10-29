"""PURPLE: Purity and ploidy estimates for somatic tumor/normal samples

https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator
"""
import csv
import os
import shutil

import toolz as tz

from bcbio import utils
from bcbio.log import logger
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
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
    purple_out = _run_purple(paired, het_file, depth_file, work_dir)
    out = []
    if paired.normal_data:
        out.append(paired.normal_data)
    if "sv" not in paired.tumor_data:
        paired.tumor_data["sv"] = []
    paired.tumor_data["sv"].append(purple_out)
    out.append(paired.tumor_data)
    return out

def _run_purple(paired, het_file, depth_file, work_dir):
    """Run PURPLE with pre-calculated AMBER and COBALT compatible inputs.

    XXX Need to add output conversion into VCF for standard formats
    """
    purple_dir = utils.safe_makedir(os.path.join(work_dir, "purple"))
    out_file = os.path.join(purple_dir, "%s.purple.cnv" % dd.get_sample_name(paired.tumor_data))
    if not utils.file_exists(out_file):
        with file_transaction(paired.tumor_data, out_file) as tx_out_file:
            cmd = ["PURPLE", "-amber", os.path.dirname(het_file), "-baf", het_file,
                   "-cobalt", os.path.dirname(depth_file),
                   "-gc_profile", dd.get_variation_resources(paired.tumor_data)["gc_profile"],
                   "-output_dir", os.path.dirname(tx_out_file),
                   "-ref_genome", "hg38" if dd.get_genome_build(paired.tumor_data) == "hg38" else "hg19",
                   "-run_dir", work_dir,
                   "-threads", dd.get_num_cores(paired.tumor_data),
                   "-tumor_sample", dd.get_sample_name(paired.tumor_data),
                   "-ref_sample", dd.get_sample_name(paired.normal_data)]
            # Avoid X11 display errors when writing plots
            cmd = "unset DISPLAY && %s" % " ".join([str(x) for x in cmd])
            do.run(cmd, "PURPLE: purity and ploidy estimation")
            for f in os.listdir(os.path.dirname(tx_out_file)):
                if f != os.path.basename(tx_out_file):
                    shutil.move(os.path.join(os.path.dirname(tx_out_file), f),
                                os.path.join(purple_dir, f))
    out_file_export = os.path.join(purple_dir, "%s-purple-cnv.tsv" % (dd.get_sample_name(paired.tumor_data)))
    if not utils.file_exists(out_file_export):
        utils.symlink_plus(out_file, out_file_export)
    out = {"variantcaller": "purple", "call_file": out_file_export,
           "plot": {}, "metrics": {}}
    for name, ext in [("copy_number", "copyNumber"), ("minor_allele", "minor_allele"), ("variant", "variant")]:
        plot_file = os.path.join(purple_dir, "plot", "%s.%s.png" % (dd.get_sample_name(paired.tumor_data), ext))
        if os.path.exists(plot_file):
            out["plot"][name] = plot_file
    purity_file = os.path.join(purple_dir, "%s.purple.purity" % dd.get_sample_name(paired.tumor_data))
    with open(purity_file) as in_handle:
        header = in_handle.readline().replace("#", "").split("\t")
        vals = in_handle.readline().split("\t")
        for h, v in zip(header, vals):
            try:
                v = float(v)
            except ValueError:
                pass
            out["metrics"][h] = v
    return out

class AmberWriter:
    def __init__(self, out_handle):
        self.writer = csv.writer(out_handle, delimiter="\t")

    def write_header(self):
        self.writer.writerow(["Chromosome", "Position", "TumorBAF", "TumorModifiedBAF", "TumorDepth",
                                "NormalBAF", "NormalModifiedBAF", "NormalDepth"])

    def _normalize_baf(self, baf):
        """Provide normalized BAF in the same manner as Amber, relative to het.

        https://github.com/hartwigmedical/hmftools/blob/637e3db1a1a995f4daefe2d0a1511a5bdadbeb05/hmf-common/src7/main/java/com/hartwig/hmftools/common/amber/AmberBAF.java#L16
        """
        if baf is None:
            baf = 0.0
        return 0.5 + abs(baf - 0.5)

    def write_row(self, rec, stats):
        if stats["normal"]["freq"] is not None and stats["normal"]["depth"] is not None:
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
                                         work_dir, paired, AmberWriter)
    amber_dir = utils.safe_makedir(os.path.join(work_dir, "amber"))
    out_file = os.path.join(amber_dir, "%s.amber.baf" % dd.get_sample_name(paired.tumor_data))
    utils.symlink_plus(prep_file, out_file)
    pcf_file = out_file + ".pcf"
    if not utils.file_exists(pcf_file):
        with file_transaction(paired.tumor_data, pcf_file) as tx_out_file:
            r_file = os.path.join(os.path.dirname(tx_out_file), "bafSegmentation.R")
            with open(r_file, "w") as out_handle:
                out_handle.write(_amber_seg_script)
            cmd = "%s && %s --no-environ %s %s %s" % (utils.get_R_exports(), utils.Rscript_cmd(), r_file,
                                                      out_file, pcf_file)
            do.run(cmd, "PURPLE: AMBER baf segmentation")
    return out_file

# BAF segmentation with copynumber from AMBER
# https://github.com/hartwigmedical/hmftools/blob/master/amber/src/main/resources/r/bafSegmentation.R
_amber_seg_script = """
# Parse the arguments
args <- commandArgs(trailing=T)
bafFile <- args[1]
pcfFile   <- args[2]

library(copynumber)
baf <- read.table(bafFile, header=TRUE)
chromosomeLevels = levels(baf$Chromosome)
chromosomePrefix = ""
if (any(grepl("chr", chromosomeLevels, ignore.case = T))) {
    chromosomePrefix = substr(chromosomeLevels[1], 1, 3)
}

baf <- baf[,c("Chromosome","Position","TumorModifiedBAF")]
baf$Chromosome <- gsub(chromosomePrefix, "", baf$Chromosome, ignore.case = T)
baf.seg<-pcf(baf,verbose=FALSE,gamma=100,kmin=1)
baf.seg$chrom = paste0(chromosomePrefix, baf.seg$chrom)
write.table(baf.seg, file = pcfFile, row.names = F, sep = "\t", quote = F)
"""

def _run_cobalt(paired, work_dir):
    """Run Cobalt for counting read depth across genomic windows.

    PURPLE requires even 1000bp windows so use integrated counting solution
    directly rather than converting from CNVkit calculations. If this approach
    is useful should be moved upstream to be available to other tools as
    an input comparison.

    https://github.com/hartwigmedical/hmftools/tree/master/count-bam-lines
    """
    cobalt_dir = utils.safe_makedir(os.path.join(work_dir, "cobalt"))
    out_file = os.path.join(cobalt_dir, "%s.cobalt" % dd.get_sample_name(paired.tumor_data))
    if not utils.file_exists(out_file):
        with file_transaction(paired.tumor_data, out_file) as tx_out_file:
            cmd = ["COBALT", "-reference", paired.normal_name, "-reference_bam", paired.normal_bam,
                   "-tumor", paired.tumor_name, "-tumor_bam", paired.tumor_bam,
                   "-threads", dd.get_num_cores(paired.tumor_data),
                   "-output_dir", os.path.dirname(tx_out_file),
                   "-gc_profile", dd.get_variation_resources(paired.tumor_data)["gc_profile"]]
            cmd = "%s && %s" % (utils.get_R_exports(), " ".join([str(x) for x in cmd]))
            do.run(cmd, "PURPLE: COBALT read depth normalization")
            for f in os.listdir(os.path.dirname(tx_out_file)):
                if f != os.path.basename(tx_out_file):
                    shutil.move(os.path.join(os.path.dirname(tx_out_file), f),
                                os.path.join(cobalt_dir, f))
    return out_file

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
