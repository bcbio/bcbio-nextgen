"""Estimate tumor purity and frequency using copy number and allele freqencies with BubbleTree.

http://www.bioconductor.org/packages/release/bioc/html/BubbleTree.html
http://www.bioconductor.org/packages/release/bioc/vignettes/BubbleTree/inst/doc/BubbleTree-vignette.html
"""
import collections
import csv
import os
import re
import subprocess

import numpy as np
import pysam
import toolz as tz

from bcbio import install, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.heterogeneity import chromhacks
from bcbio.structural import shared
from bcbio.variation import bedutils

def run(vrn_info, calls_by_name, somatic_info):
    """Run BubbleTree given variant calls, CNVs and somatic
    """
    assert "cnvkit" in calls_by_name, "BubbleTree only currently support CNVkit"
    cnv_info = calls_by_name["cnvkit"]
    work_dir = _cur_workdir(somatic_info.tumor_data)
    vcf_csv = _prep_vrn_file(vrn_info["vrn_file"], vrn_info["variantcaller"], cnv_info["cns"],
                             work_dir, calls_by_name, somatic_info)
    cnv_csv = _prep_cnv_file(cnv_info["cns"], cnv_info["variantcaller"], calls_by_name, work_dir,
                             somatic_info.tumor_data)
    return _run_bubbletree(vcf_csv, cnv_csv, somatic_info.tumor_data, somatic_info.normal_bam is not None)

def _run_bubbletree(vcf_csv, cnv_csv, data, has_normal=True):
    """Create R script and run on input data
    """
    local_sitelib = os.path.join(install.get_defaults().get("tooldir", "/usr/local"),
                                 "lib", "R", "site-library")
    base = utils.splitext_plus(vcf_csv)[0]
    r_file = "%s-run.R" % base
    bubbleplot_out = "%s-bubbleplot.pdf" % base
    trackplot_out = "%s-trackplot.pdf" % base
    calls_out = "%s-calls.rds" % base
    freqs_out = "%s-bubbletree_prevalence.txt" % base
    sample = dd.get_sample_name(data)
    # BubbleTree has some internal hardcoded paramters that assume a smaller
    # distribution of log2 scores. This is not true for tumor-only calls and
    # normal contamination, so we scale the calculations to actually get calls.
    # Need a better long term solution with flexible parameters.
    lrr_scale = 1.0 if has_normal else 10.0
    with open(r_file, "w") as out_handle:
        out_handle.write(_script.format(**locals()))
    if not utils.file_exists(freqs_out):
        try:
            do.run([utils.Rscript_cmd(), r_file], "Assess heterogeneity with BubbleTree")
        except subprocess.CalledProcessError, msg:
            if _allowed_bubbletree_errorstates(str(msg)):
                with open(freqs_out, "w") as out_handle:
                    out_handle.write('bubbletree failed:\n %s"\n' % (str(msg)))
            else:
                logger.exception()
                raise
    return {"caller": "bubbletree",
            "report": freqs_out,
            "plots": {"bubble": bubbleplot_out, "track": trackplot_out}}

def _allowed_bubbletree_errorstates(msg):
    allowed = ["Error in p[i, ] : subscript out of bounds",
               "replacement has .* rows, data has"]
    return any([len(re.findall(m, msg)) > 0 for m in allowed])

def _cns_to_coords(line):
    chrom, start, end = line.split()[:3]
    return (chrom, start, end)

def _prep_cnv_file(cns_file, svcaller, calls_by_name, work_dir, data):
    """Create a CSV file of CNV calls with log2 and number of marks.
    """
    # in_file = theta.subset_by_supported(cns_file, _cns_to_coords, calls_by_name, work_dir, data,
    #                                     headers=("chromosome", "#"))
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
                    reader.next()  # header
                    for chrom, start, end, _, log2, probes in (xs[:6] for xs in reader):
                        if chromhacks.is_autosomal(chrom):
                            writer.writerow([_to_ucsc_style(chrom), start, end, probes, log2])
    return out_file

def _prep_vrn_file(in_file, vcaller, seg_file, work_dir, calls_by_name, somatic_info):
    """Select heterozygous variants in the normal sample with sufficient depth.
    """
    data = somatic_info.tumor_data
    params = {"min_freq": 0.4,
              "max_freq": 0.6,
              "tumor_only": {"min_freq": 0.15, "max_freq": 0.85},
              "min_depth": 15,
              "hetblock": {"min_alleles": 25,
                           "allowed_misses": 2}}
    out_file = os.path.join(work_dir, "%s-%s-prep.csv" % (utils.splitext_plus(os.path.basename(in_file))[0],
                                                          vcaller))
    if not utils.file_uptodate(out_file, in_file):
        ready_bed = _identify_heterogeneity_blocks_seg(in_file, seg_file, params, work_dir, somatic_info)
        sub_file = _create_subset_file(in_file, ready_bed, work_dir, data)
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                writer = csv.writer(out_handle)
                writer.writerow(["chrom", "start", "end", "freq"])
                bcf_in = pysam.VariantFile(sub_file)
                for rec in bcf_in:
                    tumor_freq = _is_possible_loh(rec, params, somatic_info)
                    if chromhacks.is_autosomal(rec.chrom) and tumor_freq is not None:
                        writer.writerow([_to_ucsc_style(rec.chrom), rec.start, rec.stop, tumor_freq])
    return out_file

def _identify_heterogeneity_blocks_seg(in_file, seg_file, params, work_dir, somatic_info):
    """Identify heterogeneity blocks corresponding to segmentation from CNV input file.
    """
    def _segment_by_cns(target_chrom, freqs, coords):
        with open(seg_file) as in_handle:
            reader = csv.reader(in_handle, dialect="excel-tab")
            reader.next()  # header
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
                    if sum(rec.samples[somatic_info.tumor_name].allele_indices) > 0:
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

def _passes_plus_germline(rec):
    """Check if a record passes filters (but might be germline -- labelled with REJECT).
    """
    allowed = set(["PASS", "REJECT"])
    filters = [x for x in rec.filter.keys() if x not in allowed]
    return len(filters) > 0

def _is_biallelic_snp(rec):
    return _is_snp(rec) and len(rec.alts) == 1

def _is_snp(rec):
    return max([len(x) for x in rec.alleles]) == 1

def _tumor_normal_stats(rec, somatic_info):
    """Retrieve depth and frequency of tumor and normal samples.
    """
    out = {"normal": {"depth": None, "freq": None},
           "tumor": {"depth": 0, "freq": None}}
    for name, sample in rec.samples.items():
        alt, depth = sample_alt_and_depth(sample)
        if alt is not None and depth is not None and depth > 0:
            freq = float(alt) / float(depth)
            if name == somatic_info.normal_name:
                key = "normal"
            elif name == somatic_info.tumor_name:
                key = "tumor"
            out[key]["freq"] = freq
            out[key]["depth"] = depth
    return out

def _is_possible_loh(rec, params, somatic_info):
    """Check if the VCF record is a het in the normal with sufficient support.

    Only returns SNPs, since indels tend to have less precise frequency measurements.
    """
    if _is_biallelic_snp(rec) and _passes_plus_germline(rec):
        stats = _tumor_normal_stats(rec, somatic_info)
        depths = [tz.get_in([x, "depth"], stats) for x in ["normal", "tumor"]]
        depths = [d for d in depths if d is not None]
        normal_freq = tz.get_in(["normal", "freq"], stats)
        tumor_freq = tz.get_in(["tumor", "freq"], stats)
        if all([d > params["min_depth"] for d in depths]):
            if normal_freq is not None:
                if normal_freq >= params["min_freq"] and normal_freq <= params["max_freq"]:
                    return stats["tumor"]["freq"]
            elif (tumor_freq >= params["tumor_only"]["min_freq"] and
                    tumor_freq <= params["tumor_only"]["max_freq"]):
                return stats["tumor"]["freq"]

def sample_alt_and_depth(sample):
    """Flexibly get ALT allele and depth counts, handling FreeBayes, MuTect and other cases.
    """
    if "AD" in sample:
        all_counts = [int(x) for x in sample["AD"]]
        alt_counts = sum(all_counts[1:])
        depth = sum(all_counts)
    elif "AO" in sample and sample.get("RO") is not None:
        alts = sample["AO"]
        if not isinstance(alts, (list, tuple)):
            alts = [alts]
        alt_counts = sum([int(x) for x in alts])
        depth = alt_counts + int(sample["RO"])
    else:
        alt_counts = None
    if alt_counts is None:
        return None, None
    else:
        return alt_counts, depth

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
        if _is_possible_loh(rec, params, somatic(sys.argv[2], sys.argv[3])):
            print rec.filter.keys(), len(rec.filter)

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
"""
