"""PURPLE: Purity and ploidy estimates for somatic tumor/normal samples

https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator
"""
import csv
import os
import re
import shutil
import subprocess

import toolz as tz

from bcbio import broad, utils
from bcbio.log import logger
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.structural import titancna
from bcbio.variation import vcfutils

def run(items):
    paired = vcfutils.get_paired(items)
    if not paired or not paired.normal_name:
        logger.info("Skipping PURPLE; need tumor/normal somatic calls in batch: %s" %
                    " ".join([dd.get_sample_name(d) for d in items]))
        return items
    work_dir = _sv_workdir(paired.tumor_data)
    from bcbio import heterogeneity
    vrn_files = heterogeneity.get_variants(paired.tumor_data, include_germline=False)
    het_file = _amber_het_file("pon", vrn_files, work_dir, paired)
    depth_file = _run_cobalt(paired, work_dir)
    purple_out = _run_purple(paired, het_file, depth_file, vrn_files, work_dir)
    out = []
    if paired.normal_data:
        out.append(paired.normal_data)
    if "sv" not in paired.tumor_data:
        paired.tumor_data["sv"] = []
    paired.tumor_data["sv"].append(purple_out)
    out.append(paired.tumor_data)
    return out

def _get_jvm_opts(out_file, data):
    """Retrieve Java options, adjusting memory for available cores.
    """
    resources = config_utils.get_resources("purple", data["config"])
    jvm_opts = resources.get("jvm_opts", ["-Xms750m", "-Xmx3500m"])
    jvm_opts = config_utils.adjust_opts(jvm_opts, {"algorithm": {"memory_adjust":
                                                                 {"direction": "increase",
                                                                  "maximum": "30000M",
                                                                  "magnitude": dd.get_cores(data)}}})
    jvm_opts += broad.get_default_jvm_opts(os.path.dirname(out_file))
    return jvm_opts

def _run_purple(paired, het_file, depth_file, vrn_files, work_dir):
    """Run PURPLE with pre-calculated AMBER and COBALT compatible inputs.
    """
    purple_dir = utils.safe_makedir(os.path.join(work_dir, "purple"))
    out_file = os.path.join(purple_dir, "%s.purple.cnv" % dd.get_sample_name(paired.tumor_data))
    if not utils.file_exists(out_file):
        with file_transaction(paired.tumor_data, out_file) as tx_out_file:
            cmd = ["PURPLE"] + _get_jvm_opts(tx_out_file, paired.tumor_data) + \
                  ["-amber", os.path.dirname(het_file), "-baf", het_file,
                   "-cobalt", os.path.dirname(depth_file),
                   "-gc_profile", dd.get_variation_resources(paired.tumor_data)["gc_profile"],
                   "-output_dir", os.path.dirname(tx_out_file),
                   "-ref_genome", "hg38" if dd.get_genome_build(paired.tumor_data) == "hg38" else "hg19",
                   "-run_dir", work_dir,
                   "-threads", dd.get_num_cores(paired.tumor_data),
                   "-tumor_sample", dd.get_sample_name(paired.tumor_data),
                   "-ref_sample", dd.get_sample_name(paired.normal_data)]
            if vrn_files:
                cmd += ["-somatic_vcf", vrn_files[0]["vrn_file"]]
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
           "vrn_file": titancna.to_vcf(out_file_export, "PURPLE", _get_header, _export_to_vcf,
                                       paired.tumor_data),
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

def _normalize_baf(baf):
    """Provide normalized BAF in the same manner as Amber, relative to het.

    https://github.com/hartwigmedical/hmftools/blob/637e3db1a1a995f4daefe2d0a1511a5bdadbeb05/hmf-common/src7/main/java/com/hartwig/hmftools/common/amber/AmberBAF.java#L16
    """
    if baf is None:
        baf = 0.0
    return 0.5 + abs(baf - 0.5)

def _counts_to_amber(t_vals, n_vals):
    """Converts a line of CollectAllelicCounts into AMBER line.
    """
    t_depth = int(t_vals["REF_COUNT"]) + int(t_vals["ALT_COUNT"])
    n_depth = int(n_vals["REF_COUNT"]) + int(n_vals["ALT_COUNT"])
    if n_depth > 0 and t_depth > 0:
        t_baf = float(t_vals["ALT_COUNT"]) / float(t_depth)
        n_baf = float(n_vals["ALT_COUNT"]) / float(n_depth)
        return [t_vals["CONTIG"], t_vals["POSITION"], t_baf, _normalize_baf(t_baf), t_depth,
                n_baf, _normalize_baf(n_baf), n_depth]

def _count_files_to_amber(tumor_counts, normal_counts, work_dir, data):
    """Converts tumor and normal counts from GATK CollectAllelicCounts into Amber format.
    """
    amber_dir = utils.safe_makedir(os.path.join(work_dir, "amber"))
    out_file = os.path.join(amber_dir, "%s.amber.baf" % dd.get_sample_name(data))

    if not utils.file_uptodate(out_file, tumor_counts):
        with file_transaction(data, out_file) as tx_out_file:
            with open(tumor_counts) as tumor_handle:
                with open(normal_counts) as normal_handle:
                    with open(tx_out_file, "w") as out_handle:
                        writer = csv.writer(out_handle, delimiter="\t")
                        writer.writerow(["Chromosome", "Position", "TumorBAF", "TumorModifiedBAF", "TumorDepth",
                                         "NormalBAF", "NormalModifiedBAF", "NormalDepth"])
                        header = None
                        for t, n in zip(tumor_handle, normal_handle):
                            if header is None and t.startswith("CONTIG"):
                                header = t.strip().split()
                            elif header is not None:
                                t_vals = dict(zip(header, t.strip().split()))
                                n_vals = dict(zip(header, n.strip().split()))
                                amber_line = _counts_to_amber(t_vals, n_vals)
                                if amber_line:
                                    writer.writerow(amber_line)
    return out_file

class AmberWriter:
    def __init__(self, out_handle):
        self.writer = csv.writer(out_handle, delimiter="\t")

    def write_header(self):
        self.writer.writerow(["Chromosome", "Position", "TumorBAF", "TumorModifiedBAF", "TumorDepth",
                                "NormalBAF", "NormalModifiedBAF", "NormalDepth"])

    def write_row(self, rec, stats):
        if stats["normal"]["freq"] is not None and stats["normal"]["depth"] is not None:
            self.writer.writerow([rec.chrom, rec.pos,
                                  stats["tumor"]["freq"], _normalize_baf(stats["tumor"]["freq"]),
                                  stats["tumor"]["depth"],
                                  stats["normal"]["freq"], _normalize_baf(stats["normal"]["freq"]),
                                  stats["normal"]["depth"]])

def _amber_het_file(method, vrn_files, work_dir, paired):
    """Create file of BAFs in normal heterozygous positions compatible with AMBER.

    Two available methods:
      - pon -- Use panel of normals with likely heterozygous sites.
      - variants -- Use pre-existing variant calls, filtered to likely heterozygotes.

    https://github.com/hartwigmedical/hmftools/tree/master/amber
    https://github.com/hartwigmedical/hmftools/blob/637e3db1a1a995f4daefe2d0a1511a5bdadbeb05/hmf-common/src/test/resources/amber/new.amber.baf
    """
    assert vrn_files, "Did not find compatible variant calling files for PURPLE inputs"
    from bcbio.heterogeneity import bubbletree

    if method == "variants":
        amber_dir = utils.safe_makedir(os.path.join(work_dir, "amber"))
        out_file = os.path.join(amber_dir, "%s.amber.baf" % dd.get_sample_name(paired.tumor_data))
        prep_file = bubbletree.prep_vrn_file(vrn_files[0]["vrn_file"], vrn_files[0]["variantcaller"],
                                             work_dir, paired, AmberWriter)
        utils.symlink_plus(prep_file, out_file)
        pcf_file = out_file + ".pcf"
        if not utils.file_exists(pcf_file):
            with file_transaction(paired.tumor_data, pcf_file) as tx_out_file:
                r_file = os.path.join(os.path.dirname(tx_out_file), "bafSegmentation.R")
                with open(r_file, "w") as out_handle:
                    out_handle.write(_amber_seg_script)
                cmd = "%s && %s --vanilla %s %s %s" % (utils.get_R_exports(), utils.Rscript_cmd(), r_file,
                                                          out_file, pcf_file)
                do.run(cmd, "PURPLE: AMBER baf segmentation")
    else:
        assert method == "pon"
        out_file = _run_amber(paired, work_dir)
    return out_file

def _run_amber(paired, work_dir, lenient=False):
    """AMBER: calculate allele frequencies at likely heterozygous sites.

    lenient flag allows amber runs on small test sets.
    """
    amber_dir = utils.safe_makedir(os.path.join(work_dir, "amber"))
    out_file = os.path.join(amber_dir, "%s.amber.baf" % dd.get_sample_name(paired.tumor_data))
    if not utils.file_exists(out_file) or not utils.file_exists(out_file + ".pcf"):
        with file_transaction(paired.tumor_data, out_file) as tx_out_file:
            key = "germline_het_pon"
            het_bed = tz.get_in(["genome_resources", "variation", key], paired.tumor_data)
            cmd = ["AMBER"] + _get_jvm_opts(tx_out_file, paired.tumor_data) + \
                  ["-threads", dd.get_num_cores(paired.tumor_data),
                   "-tumor", dd.get_sample_name(paired.tumor_data),
                   "-tumor_bam", dd.get_align_bam(paired.tumor_data),
                   "-reference", dd.get_sample_name(paired.normal_data),
                   "-reference_bam", dd.get_align_bam(paired.normal_data),
                   "-ref_genome", dd.get_ref_file(paired.tumor_data),
                   "-bed", het_bed,
                   "-output_dir", os.path.dirname(tx_out_file)]
            if lenient:
                cmd += ["-max_het_af_percent", "1.0"]
            try:
                do.run(cmd, "PURPLE: AMBER baf generation")
            except subprocess.CalledProcessError as msg:
                if not lenient and _amber_allowed_errors(str(msg)):
                    return _run_amber(paired, work_dir, True)
            for f in os.listdir(os.path.dirname(tx_out_file)):
                if f != os.path.basename(tx_out_file):
                    shutil.move(os.path.join(os.path.dirname(tx_out_file), f),
                                os.path.join(amber_dir, f))
    return out_file

def _amber_allowed_errors(msg):
    allowed = ["R execution failed. Unable to complete segmentation."]
    return any([len(re.findall(m, msg)) > 0 for m in allowed])

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
            cmd = ["COBALT"] + _get_jvm_opts(tx_out_file, paired.tumor_data) + \
                  ["-reference", paired.normal_name, "-reference_bam", paired.normal_bam,
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
        raise NotImplementedError
    return out_file

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(dd.get_work_dir(data), "structural",
                                           dd.get_sample_name(data), "purple"))

# ## VCF output

def _get_header(in_handle):
    return in_handle.readline().replace("#", "").strip().split(), in_handle

def _export_to_vcf(cur):
    """Convert PURPLE custom output into VCF.
    """
    if float(cur["copyNumber"]) > 2.0:
        svtype = "DUP"
    elif float(cur["copyNumber"]) < 2.0:
        svtype = "DEL"
    else:
        svtype = None
    if svtype:
        info = ["END=%s" % cur["end"], "SVLEN=%s" % (int(cur["end"]) - int(cur["start"])),
                "SVTYPE=%s" % svtype, "CN=%s" % cur["copyNumber"], "PROBES=%s" % cur["depthWindowCount"]]
        return [cur["chromosome"], cur["start"], ".", "N", "<%s>" % svtype, ".", ".",
                ";".join(info), "GT", "0/1"]
