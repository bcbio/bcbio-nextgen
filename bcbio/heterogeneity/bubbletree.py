"""Estimate tumor purity and frequency using copy number and allele freqencies with BubbleTree.

http://www.bioconductor.org/packages/release/bioc/html/BubbleTree.html
http://www.bioconductor.org/packages/release/bioc/vignettes/BubbleTree/inst/doc/BubbleTree-vignette.html
"""
import csv
import os

from pysam import VariantFile

from bcbio import install, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.structural import regions

def run(vrn_info, cnvs_by_name, somatic_info):
    """Run BubbleTree given variant calls, CNVs and somatic
    """
    work_dir = _cur_workdir(somatic_info.tumor_data)
    vcf_csv = _prep_vrn_file(vrn_info["vrn_file"], vrn_info["variantcaller"], work_dir, somatic_info)
    assert "cnvkit" in cnvs_by_name, "BubbleTree only currently support CNVkit"
    cnv_info = cnvs_by_name["cnvkit"]
    cnv_csv = _prep_cnv_file(cnv_info["cns"], cnv_info["variantcaller"], work_dir, somatic_info.tumor_data)
    _run_bubbletree(vcf_csv, cnv_csv, somatic_info.tumor_data)

def _run_bubbletree(vcf_csv, cnv_csv, data):
    """Create R script and run on input data
    """
    local_sitelib = os.path.join(install.get_defaults().get("tooldir", "/usr/local"),
                                 "lib", "R", "site-library")
    base = utils.splitext_plus(vcf_csv)[0]
    r_file = "%s-run.R" % base
    bubbles_out = "%s-bubbles.pdf" % base
    prev_model_out = "%s-bubbletree_prev_model.pdf" % base
    freqs_out = "%s-bubbletree_prevalence.txt" % base
    with open(r_file, "w") as out_handle:
        out_handle.write(_script.format(**locals()))
    if not utils.file_exists(freqs_out):
        do.run(["Rscript", r_file], "Assess heterogeneity with BubbleTree")

def _prep_cnv_file(in_file, svcaller, work_dir, data):
    """Create a CSV file of CNV calls with log2 and number of marks.
    """
    out_file = os.path.join(work_dir, "%s-%s-prep.csv" % (utils.splitext_plus(os.path.basename(in_file))[0],
                                                          svcaller))
    autosomal_chroms = _get_autosomal_chroms()
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(in_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    reader = csv.reader(in_handle, dialect="excel-tab")
                    writer = csv.writer(out_handle)
                    writer.writerow(["chrom", "start", "end", "num.mark", "seg.mean"])
                    reader.next()  # header
                    for chrom, start, end, _, log2, probes in reader:
                        if chrom in autosomal_chroms:
                            writer.writerow([_to_ucsc_style(chrom), start, end, probes, log2])
    return out_file

def _prep_vrn_file(in_file, vcaller, work_dir, somatic_info):
    """Select heterozygous variants in the normal sample with sufficient depth.
    """
    data = somatic_info.tumor_data
    params = {"min_freq": 0.4,
              "max_freq": 0.6,
              "min_depth": 15}
    autosomal_chroms = _get_autosomal_chroms()
    out_file = os.path.join(work_dir, "%s-%s-prep.csv" % (utils.splitext_plus(os.path.basename(in_file))[0],
                                                          vcaller))
    if not utils.file_uptodate(out_file, in_file):
        sub_file = _create_subset_file(in_file, work_dir, data)
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                writer = csv.writer(out_handle)
                writer.writerow(["chrom", "start", "end", "freq"])
                bcf_in = VariantFile(sub_file)
                for rec in bcf_in:
                    tumor_freq = _is_possible_loh(rec, params, somatic_info)
                    if rec.chrom in autosomal_chroms and tumor_freq is not None:
                        writer.writerow([_to_ucsc_style(rec.chrom), rec.start, rec.stop, tumor_freq])
    return out_file

def _create_subset_file(in_file, work_dir, data):
    """Subset the VCF to a set of smaller regions, matching what was used for CNV calling.
    """
    out_file = os.path.join(work_dir, "%s-orig.bcf" % utils.splitext_plus(os.path.basename(in_file))[0])
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(data, out_file) as tx_out_file:
            region_bed = regions.get_sv_bed(data)
            if not region_bed:
                region_bed = regions.get_sv_bed(data, "transcripts1e4", work_dir)
            cmd = "bcftools view -R {region_bed} -o {tx_out_file} -O b {in_file}"
            do.run(cmd.format(**locals()), "Extract SV only regions for BubbleTree")
    return out_file

def _to_ucsc_style(chrom):
    """BubbleTree assumes hg19 UCSC style chromosome inputs.
    """
    return "chr%s" % chrom if not str(chrom).startswith("chr") else chrom

def _get_autosomal_chroms():
    """Hack to only use autosomal chromosomes that should generalize to any species with numeric chroms.
    """
    max_chrom = 1000
    return set([str(x) for x in range(1, max_chrom)] + ["chr%s" % x for x in range(1, max_chrom)])

def _is_snp(rec):
    return max([len(x) for x in rec.alleles]) == 1

def _is_possible_loh(rec, params, somatic_info):
    """Check if the VCF record is a het in the normal with sufficient support.

    Only returns SNPs, since indels tend to have less precise frequency measurements.
    """
    normal_good, tumor_good = False, False
    if _is_snp(rec):
        for name, sample in rec.samples.items():
            alt, depth = sample_alt_and_depth(sample)
            if alt is not None and depth is not None and depth > 0:
                freq = float(alt) / float(depth)
                if name == somatic_info.normal_name:
                    normal_good = (depth >= params["min_depth"] and
                                   (freq >= params["min_freq"] and freq <= params["max_freq"]))
                elif name == somatic_info.tumor_name:
                    tumor_good = (depth >= params["min_depth"] and
                                  (freq < params["min_freq"] or freq > params["max_freq"]))
                    tumor_freq = freq
    if normal_good and tumor_good:
        return tumor_freq

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
    import collections
    bcf_in = VariantFile(sys.argv[1])
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

vc.df = read.csv("{vcf_csv}", header=T)
vc.gr = GRanges(vc.df$chrom, IRanges(vc.df$start, vc.df$end), freq=vc.df$freq)

cnv.df = read.csv("{cnv_csv}", header=T)
cnv.gr = GRanges(cnv.df$chrom, IRanges(cnv.df$start, cnv.df$end),
                 num.mark=cnv.df$num.mark, seg.mean=cnv.df$seg.mean)
print(vc.gr)
print(cnv.gr)

rbd = getRBD(snp.gr=vc.gr, cnv.gr=cnv.gr)
pdf(file="{bubbles_out}", width=8, height=6)
plotBubbles(rbd)
dev.off()
pur = calc.prev(rbd, heurx=FALSE, modex=5, plotx="{prev_model_out}")

sink("{freqs_out}")
# tumor subclone freqencies
print(pur[[1]]$ploidy_prev)
# tumor purity
print(pur[[2]][nrow(pur[[2]]),2])
sink()
"""
