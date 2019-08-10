"""Estimate tumor purity and frequency using copy number and allele freqencies with BubbleTree.

http://www.bioconductor.org/packages/release/bioc/html/BubbleTree.html
http://www.bioconductor.org/packages/release/bioc/vignettes/BubbleTree/inst/doc/BubbleTree-vignette.html
"""
from __future__ import print_function
import collections
import csv
import os
import re
import subprocess

import numpy as np
import pysam
import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.heterogeneity import chromhacks
from bcbio.structural import shared
from bcbio.variation import bedutils

population_keys = ['AC_AFR', 'AC_AMR', 'AC_EAS', 'AC_FIN', 'AC_NFE', 'AC_OTH', 'AC_SAS']
PARAMS = {"min_freq": 0.2,
          "max_freq": 0.8,
          "tumor_only": {"min_freq": 0.10, "max_freq": 0.90},
          "min_depth": 15,
          "hetblock": {"min_alleles": 25,
                       "allowed_misses": 2}}

def run(vrn_info, calls_by_name, somatic_info, do_plots=True, handle_failures=True):
    """Run BubbleTree given variant calls, CNVs and somatic
    """
    if "seq2c" in calls_by_name:
        cnv_info = calls_by_name["seq2c"]
    elif "cnvkit" in calls_by_name:
        cnv_info = calls_by_name["cnvkit"]
    else:
        raise ValueError("BubbleTree only currently support CNVkit and Seq2c: %s" % ", ".join(calls_by_name.keys()))
    work_dir = _cur_workdir(somatic_info.tumor_data)
    class OutWriter:
        def __init__(self, out_handle):
            self.writer = csv.writer(out_handle)
        def write_header(self):
            self.writer.writerow(["chrom", "start", "end", "freq"])
        def write_row(self, rec, stats):
            self.writer.writerow([_to_ucsc_style(rec.chrom), rec.start, rec.stop, stats["tumor"]["freq"]])
    vcf_csv = prep_vrn_file(vrn_info["vrn_file"], vrn_info["variantcaller"],
                            work_dir, somatic_info, OutWriter, cnv_info["cns"])
    cnv_csv = _prep_cnv_file(cnv_info["cns"], cnv_info["variantcaller"], work_dir,
                             somatic_info.tumor_data)
    wide_lrr = cnv_info["variantcaller"] == "cnvkit" and somatic_info.normal_bam is None
    return _run_bubbletree(vcf_csv, cnv_csv, somatic_info.tumor_data, wide_lrr, do_plots,
                           handle_failures)

def _run_bubbletree(vcf_csv, cnv_csv, data, wide_lrr=False, do_plots=True,
                    handle_failures=True):
    """Create R script and run on input data

    BubbleTree has some internal hardcoded paramters that assume a smaller
    distribution of log2 scores. This is not true for tumor-only calls, so if
    we specify wide_lrr we scale the calculations to actually get calls. Need a
    better long term solution with flexible parameters.
    """
    lrr_scale = 10.0 if wide_lrr else 1.0
    local_sitelib = utils.R_sitelib()
    base = utils.splitext_plus(vcf_csv)[0]
    r_file = "%s-run.R" % base
    bubbleplot_out = "%s-bubbleplot.pdf" % base
    trackplot_out = "%s-trackplot.pdf" % base
    calls_out = "%s-calls.rds" % base
    freqs_out = "%s-bubbletree_prevalence.txt" % base
    sample = dd.get_sample_name(data)
    do_plots = "yes" if do_plots else "no"
    with open(r_file, "w") as out_handle:
        out_handle.write(_script.format(**locals()))
    if not utils.file_exists(freqs_out):
        cmd = "%s && %s --vanilla %s" % (utils.get_R_exports(), utils.Rscript_cmd(), r_file)
        try:
            do.run(cmd, "Assess heterogeneity with BubbleTree")
        except subprocess.CalledProcessError as msg:
            if handle_failures and _allowed_bubbletree_errorstates(str(msg)):
                with open(freqs_out, "w") as out_handle:
                    out_handle.write('bubbletree failed:\n %s"\n' % (str(msg)))
            else:
                logger.exception()
                raise
    return {"caller": "bubbletree",
            "report": freqs_out,
            "plot": {"bubble": bubbleplot_out, "track": trackplot_out}}

def _allowed_bubbletree_errorstates(msg):
    allowed = ["Error in p[i, ] : subscript out of bounds",
               "replacement has .* rows, data has"]
    return any([len(re.findall(m, msg)) > 0 for m in allowed])

def _cns_to_coords(line):
    chrom, start, end = line.split()[:3]
    return (chrom, start, end)

def _prep_cnv_file(cns_file, svcaller, work_dir, data):
    """Create a CSV file of CNV calls with log2 and number of marks.
    """
    in_file = cns_file
    out_file = os.path.join(work_dir, "%s-%s-prep.csv" % (utils.splitext_plus(os.path.basename(in_file))[0],
                                                          svcaller))
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(in_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    reader = csv.reader(in_handle, dialect="excel-tab")
                    writer = csv.writer(out_handle)
                    writer.writerow(["chrom", "start", "end", "num.mark", "seg.mean"])
                    header = next(reader)
                    for line in reader:
                        cur = dict(zip(header, line))
                        if chromhacks.is_autosomal(cur["chromosome"]):
                            writer.writerow([_to_ucsc_style(cur["chromosome"]), cur["start"],
                                             cur["end"], cur["probes"], cur["log2"]])
    return out_file

def prep_vrn_file(in_file, vcaller, work_dir, somatic_info, writer_class, seg_file=None, params=None):
    """Select heterozygous variants in the normal sample with sufficient depth.

    writer_class implements write_header and write_row to write VCF outputs
    from a record and extracted tumor/normal statistics.
    """
    data = somatic_info.tumor_data
    if not params:
        params = PARAMS
    out_file = os.path.join(work_dir, "%s-%s-prep.csv" % (utils.splitext_plus(os.path.basename(in_file))[0],
                                                          vcaller))
    if not utils.file_uptodate(out_file, in_file):
        # ready_bed = _identify_heterogeneity_blocks_seg(in_file, seg_file, params, work_dir, somatic_info)
        ready_bed = None
        if ready_bed and utils.file_exists(ready_bed):
            sub_file = _create_subset_file(in_file, ready_bed, work_dir, data)
        else:
            sub_file = in_file
        max_depth = max_normal_germline_depth(sub_file, params, somatic_info)
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                writer = writer_class(out_handle)
                writer.write_header()
                bcf_in = pysam.VariantFile(sub_file)
                for rec in bcf_in:
                    stats = _is_possible_loh(rec, bcf_in, params, somatic_info, max_normal_depth=max_depth)
                    if chromhacks.is_autosomal(rec.chrom) and stats is not None:
                        writer.write_row(rec, stats)
    return out_file


# thresholds for filtering normal samples based on depth
# matches those used in PURPLE AMBER caller
# https://github.com/hartwigmedical/hmftools/blob/a8c5dd2487c8c294c457eb8961c78e08c61a604a/amber/src/main/java/com/hartwig/hmftools/amber/AmberApplication.java#L41
NORMAL_FILTER_PARAMS = {"min_depth_percent": 0.5, "max_depth_percent": 1.5,
                        "min_freq_narrow": 0.4, "max_freq_narrow": 0.65}

def max_normal_germline_depth(in_file, params, somatic_info):
    """Calculate threshold for excluding potential heterozygotes based on normal depth.
    """
    bcf_in = pysam.VariantFile(in_file)
    depths = []
    for rec in bcf_in:
        stats = _is_possible_loh(rec, bcf_in, params, somatic_info)
        if tz.get_in(["normal", "depth"], stats):
            depths.append(tz.get_in(["normal", "depth"], stats))
    if depths:
        return np.median(depths) * NORMAL_FILTER_PARAMS["max_depth_percent"]

def _identify_heterogeneity_blocks_seg(in_file, seg_file, params, work_dir, somatic_info):
    """Identify heterogeneity blocks corresponding to segmentation from CNV input file.
    """
    def _segment_by_cns(target_chrom, freqs, coords):
        with open(seg_file) as in_handle:
            reader = csv.reader(in_handle, dialect="excel-tab")
            next(reader)  # header
            for cur_chrom, start, end in (xs[:3] for xs in reader):
                if cur_chrom == target_chrom:
                    block_freqs = []
                    for i, (freq, coord) in enumerate(zip(freqs, coords)):
                        if coord >= int(start) and coord < int(end):
                            block_freqs.append(freq)
                        elif coord >= int(end):
                            break
                    coords = coords[max(0, i - 1):]
                    freqs = freqs[max(0, i - 1):]
                    if len(block_freqs) > params["hetblock"]["min_alleles"]:
                        yield start, end
    return _identify_heterogeneity_blocks_shared(in_file, _segment_by_cns, params, work_dir, somatic_info)

def _identify_heterogeneity_blocks_hmm(in_file, params, work_dir, somatic_info):
    """Use a HMM to identify blocks of heterogeneity to use for calculating allele frequencies.

    The goal is to subset the genome to a more reasonable section that contains potential
    loss of heterogeneity or other allele frequency adjustment based on selection.
    """
    def _segment_by_hmm(chrom, freqs, coords):
        cur_coords = []
        for j, state in enumerate(_predict_states(freqs)):
            if state == 0:  # heterozygote region
                if len(cur_coords) == 0:
                    num_misses = 0
                cur_coords.append(coords[j])
            else:
                num_misses += 1
            if num_misses > params["hetblock"]["allowed_misses"]:
                if len(cur_coords) >= params["hetblock"]["min_alleles"]:
                    yield min(cur_coords), max(cur_coords)
                cur_coords = []
        if len(cur_coords) >= params["hetblock"]["min_alleles"]:
            yield min(cur_coords), max(cur_coords)
    return _identify_heterogeneity_blocks_shared(in_file, _segment_by_hmm, params, work_dir, somatic_info)

def _identify_heterogeneity_blocks_shared(in_file, segment_fn, params, work_dir, somatic_info):
    """Identify heterogeneity blocks corresponding to segmentation from CNV input file.
    """
    out_file = os.path.join(work_dir, "%s-hetblocks.bed" % utils.splitext_plus(os.path.basename(in_file))[0])
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(somatic_info.tumor_data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                for chrom, freqs, coords in _freqs_by_chromosome(in_file, params, somatic_info):
                    for start, end in segment_fn(chrom, freqs, coords):
                        out_handle.write("%s\t%s\t%s\n" % (chrom, start, end))
    return out_file

def _predict_states(freqs):
    """Use frequencies to predict states across a chromosome.

    Normalize so heterozygote blocks are assigned state 0 and homozygous
    are assigned state 1.
    """
    from hmmlearn import hmm
    freqs = np.column_stack([np.array(freqs)])
    model = hmm.GaussianHMM(2, covariance_type="full")
    model.fit(freqs)
    states = model.predict(freqs)
    freqs_by_state = collections.defaultdict(list)
    for i, state in enumerate(states):
        freqs_by_state[state].append(freqs[i])
    if np.median(freqs_by_state[0]) > np.median(freqs_by_state[1]):
        states = [0 if s == 1 else 1 for s in states]
    return states

def _freqs_by_chromosome(in_file, params, somatic_info):
    """Retrieve frequencies across each chromosome as inputs to HMM.
    """
    freqs = []
    coords = []
    cur_chrom = None
    with pysam.VariantFile(in_file) as bcf_in:
        for rec in bcf_in:
            if _is_biallelic_snp(rec) and _passes_plus_germline(rec) and chromhacks.is_autosomal(rec.chrom):
                if cur_chrom is None or rec.chrom != cur_chrom:
                    if cur_chrom and len(freqs) > 0:
                        yield cur_chrom, freqs, coords
                    cur_chrom = rec.chrom
                    freqs = []
                    coords = []
                stats = _tumor_normal_stats(rec, somatic_info)
                if tz.get_in(["tumor", "depth"], stats, 0) > params["min_depth"]:
                    # not a ref only call
                    if len(rec.samples) == 0 or sum(rec.samples[somatic_info.tumor_name].allele_indices) > 0:
                        freqs.append(tz.get_in(["tumor", "freq"], stats))
                        coords.append(rec.start)
        if cur_chrom and len(freqs) > 0:
            yield cur_chrom, freqs, coords

def _create_subset_file(in_file, het_region_bed, work_dir, data):
    """Subset the VCF to a set of pre-calculated smaller regions.
    """
    cnv_regions = shared.get_base_cnv_regions(data, work_dir)
    region_bed = bedutils.intersect_two(het_region_bed, cnv_regions, work_dir, data)
    out_file = os.path.join(work_dir, "%s-origsubset.bcf" % utils.splitext_plus(os.path.basename(in_file))[0])
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(data, out_file) as tx_out_file:
            regions = ("-R %s" % region_bed) if utils.file_exists(region_bed) else ""
            cmd = "bcftools view {regions} -o {tx_out_file} -O b {in_file}"
            do.run(cmd.format(**locals()), "Extract regions for BubbleTree frequency determination")
    return out_file

def _to_ucsc_style(chrom):
    """BubbleTree assumes hg19 UCSC style chromosome inputs.
    """
    return "chr%s" % chrom if not str(chrom).startswith("chr") else chrom

def is_info_germline(rec):
    """Check if a variant record is germline based on INFO attributes.

    Works with VarDict's annotation of STATUS.
    """
    if hasattr(rec, "INFO"):
        status = rec.INFO.get("STATUS", "").lower()
    else:
        status = rec.info.get("STATUS", "").lower()
    return status == "germline" or status.find("loh") >= 0

def _passes_plus_germline(rec, use_status=False):
    """Check if a record passes filters (but might be germline -- labelled with REJECT).
    """
    if use_status and is_info_germline(rec):
        return True
    allowed = set(["PASS", "REJECT", "."])
    if hasattr(rec, "FILTER"):
        if not rec.FILTER:
            filters = []
        else:
            filters = [x for x in rec.FILTER.split(";") if x not in allowed]
    else:
        filters = [x for x in rec.filter.keys() if x not in allowed]
    return len(filters) == 0

def _is_biallelic_snp(rec):
    if hasattr(rec, "ALT"):
        return _is_snp(rec) and len(rec.ALT) == 1
    else:
        return _is_snp(rec) and len(rec.alts) == 1

def _is_snp(rec):
    if hasattr(rec, "ALT"):
        return max([len(x) for x in rec.ALT + [rec.REF]]) == 1
    else:
        return max([len(x) for x in rec.alleles]) == 1

def _tumor_normal_stats(rec, somatic_info, vcf_rec):
    """Retrieve depth and frequency of tumor and normal samples.
    """
    out = {"normal": {"alt": None, "depth": None, "freq": None},
           "tumor": {"alt": 0, "depth": 0, "freq": None}}
    if hasattr(vcf_rec, "samples"):
        samples = [(s, {}) for s in vcf_rec.samples]
        for fkey in ["AD", "AO", "RO", "AF", "DP"]:
            try:
                for i, v in enumerate(rec.format(fkey)):
                    samples[i][1][fkey] = v
            except KeyError:
                pass
    # Handle INFO only inputs
    elif len(rec.samples) == 0:
        samples = [(somatic_info.tumor_name, None)]
    else:
        samples = rec.samples.items()
    for name, sample in samples:
        alt, depth, freq = sample_alt_and_depth(rec, sample)
        if depth is not None and freq is not None:
            if name == somatic_info.normal_name:
                key = "normal"
            elif name == somatic_info.tumor_name:
                key = "tumor"
            out[key]["freq"] = freq
            out[key]["depth"] = depth
            out[key]["alt"] = alt
    return out

def _is_possible_loh(rec, vcf_rec, params, somatic_info, use_status=False, max_normal_depth=None):
    """Check if the VCF record is a het in the normal with sufficient support.

    Only returns SNPs, since indels tend to have less precise frequency measurements.
    """
    if _is_biallelic_snp(rec) and _passes_plus_germline(rec, use_status=use_status):
        stats = _tumor_normal_stats(rec, somatic_info, vcf_rec)
        depths = [tz.get_in([x, "depth"], stats) for x in ["normal", "tumor"]]
        depths = [d for d in depths if d is not None]
        normal_freq = tz.get_in(["normal", "freq"], stats)
        tumor_freq = tz.get_in(["tumor", "freq"], stats)
        if all([d > params["min_depth"] for d in depths]):
            if max_normal_depth and tz.get_in(["normal", "depth"], stats, 0) > max_normal_depth:
                return None
            if normal_freq is not None:
                if normal_freq >= params["min_freq"] and normal_freq <= params["max_freq"]:
                    return stats
            elif (tumor_freq >= params["tumor_only"]["min_freq"] and
                    tumor_freq <= params["tumor_only"]["max_freq"]):
                if (vcf_rec and not _has_population_germline(vcf_rec)) or is_population_germline(rec):
                    return stats

def _has_population_germline(rec):
    """Check if header defines population annotated germline samples for tumor only.
    """
    for k in population_keys:
        if k in rec.header.info:
            return True
    return False

def is_population_germline(rec):
    """Identify a germline calls based on annoations with ExAC or other population databases.
    """
    min_count = 50
    for k in population_keys:
        if k in rec.info:
            val = rec.info.get(k)
            if "," in val:
                val = val.split(",")[0]
            if isinstance(val, (list, tuple)):
                val = max(val)
            if int(val) > min_count:
                return True
    return False

def sample_alt_and_depth(rec, sample):
    """Flexibly get ALT allele and depth counts, handling FreeBayes, MuTect and other cases.
    """
    if sample and "AD" in sample:
        all_counts = [int(x) for x in sample["AD"]]
        alt_counts = sum(all_counts[1:])
        depth = sum(all_counts)
    elif sample and "AO" in sample and sample.get("RO") is not None:
        alts = sample["AO"]
        if not isinstance(alts, (list, tuple)):
            alts = [alts]
        alt_counts = sum([int(x) for x in alts])
        depth = alt_counts + int(sample["RO"])
    elif "DP" in rec.info and "AF" in rec.info:
        af = rec.info["AF"][0] if isinstance(rec.info["AF"], (tuple, list)) else rec.info["AF"]
        return None, rec.info["DP"], af
    else:
        alt_counts = None
    if alt_counts is None or depth is None or depth == 0:
        return None, None, None
    else:
        freq = float(alt_counts) / float(depth)
        return alt_counts, depth, freq

def _cur_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "heterogeneity",
                                           dd.get_sample_name(data), "bubbletree"))

if __name__ == "__main__":
    import sys
    bcf_in = pysam.VariantFile(sys.argv[1])
    somatic = collections.namedtuple("Somatic", "normal_name,tumor_name")
    params = {"min_freq": 0.4,
              "max_freq": 0.6,
              "min_depth": 15}
    for rec in bcf_in:
        if _is_possible_loh(rec, bcf_in, params, somatic(sys.argv[2], sys.argv[3])):
            print(rec.filter.keys(), len(rec.filter))

_script = """
.libPaths(c("{local_sitelib}"))
library(BubbleTree)
library(GenomicRanges)
library(ggplot2)

vc.df = read.csv("{vcf_csv}", header=T)
vc.gr = GRanges(vc.df$chrom, IRanges(vc.df$start, vc.df$end),
                freq=vc.df$freq, score=vc.df$freq)

cnv.df = read.csv("{cnv_csv}", header=T)
cnv.gr = GRanges(cnv.df$chrom, IRanges(cnv.df$start, cnv.df$end),
                 num.mark=cnv.df$num.mark, seg.mean=cnv.df$seg.mean,
                 score=cnv.df$seg.mean)
print(vc.gr)
print(cnv.gr)

r <- new("RBD")
rbd <- makeRBD(r, vc.gr, cnv.gr)
rbd$lrr <- rbd$lrr / {lrr_scale}
print(head(rbd))
calls <- new("BTreePredictor", rbd=rbd)
calls <- btpredict(calls)

saveRDS(calls, "{calls_out}")

purity <- calls@result$prev[1]
adj <- calls@result$ploidy.adj["adj"]
# when purity is low the calculation result is not reliable
ploidy <- (2*adj -2)/purity + 2
out <- data.frame(sample="{sample}",
                  purity=round(purity,3),
                  prevalences=paste(round(calls@result$prev,3), collapse=";"),
                  tumor_ploidy=round(ploidy,1))
write.csv(out, file="{freqs_out}", row.names=FALSE)

title <- sprintf("{sample} (%s)", info(calls))

# XXX Needs to be generalized for non-build 37/hg19 plots
# hg19.seqinfo is hardcoded in TrackPlotter but we might
# be able to work around with just an external centromere.dat import
if ("{do_plots}" == "yes") {{
  load(system.file("data", "centromere.dat.rda", package="BubbleTree"))
  load(system.file("data", "hg19.seqinfo.rda", package="BubbleTree"))
  trackplotter <- new("TrackPlotter")
  z1 <- heteroLociTrack(trackplotter,  calls@result, centromere.dat, vc.gr) + ggplot2::labs(title=title)
  z2 <- RscoreTrack(trackplotter, calls@result, centromere.dat, cnv.gr)
  t2 <- getTracks(z1, z2)
  pdf(file="{trackplot_out}", width=8, height=6)
  g <- gridExtra::grid.arrange(t2, ncol=1)
  print(g)
  dev.off()

  pdf(file="{bubbleplot_out}", width=8, height=6)
  btreeplotter <- new("BTreePlotter")
  g <- drawBTree(btreeplotter, calls@rbd.adj) + ggplot2::labs(title=title)
  print(g)
  dev.off()
}}
"""
