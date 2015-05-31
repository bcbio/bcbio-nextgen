"""Estimate tumor purity and frequency using copy number and allele freqencies with BubbleTree.

http://www.bioconductor.org/packages/release/bioc/html/BubbleTree.html
http://www.bioconductor.org/packages/release/bioc/vignettes/BubbleTree/inst/doc/BubbleTree-vignette.html
"""
import os

from pysam import VariantFile

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import regions

def run(vrn_info, cnv_info, somatic_info):
    """Run BubbleTree given variant calls, CNVs and somatic
    """
    work_dir = _cur_workdir(somatic_info.tumor_data)
    input_vcf = _prep_vrn_file(vrn_info["vrn_file"], vrn_info["variantcaller"], work_dir, somatic_info)
    print input_vcf
    print cnv_info

def _prep_vrn_file(in_file, vcaller, work_dir, somatic_info):
    """Select heterozygous variants in the normal sample with sufficient depth.
    """
    data = somatic_info.tumor_data
    params = {"min_freq": 0.4,
              "max_freq": 0.6,
              "min_depth": 15}
    out_file = os.path.join(work_dir, "%s-%s-prep.vcf" % (utils.splitext_plus(os.path.basename(in_file))[0],
                                                          vcaller))
    if not utils.file_uptodate(out_file, in_file):
        sub_file = _create_subset_file(in_file, work_dir, data)
        with file_transaction(data, out_file) as tx_out_file:
            bcf_in = VariantFile(sub_file)
            bcf_out = VariantFile(tx_out_file, "w", header=bcf_in.header)
            for rec in bcf_in:
                if _is_possible_loh(rec, params, somatic_info):
                    bcf_out.write(rec)
    return out_file

def _create_subset_file(in_file, work_dir, data):
    """Subset the VCF to a set of smaller regions, matching what was used for CNV calling.
    """
    out_file = "%s-orig.bcf" % os.path.splitext(in_file)[0]
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(data, out_file) as tx_out_file:
            region_bed = regions.get_sv_bed(data)
            if not region_bed:
                region_bed = regions.get_sv_bed(data, "transcripts1e4", work_dir)
            cmd = "vcftools view -R {region_bed} -o {tx_out_file} -O b {in_file}"
            do.run(cmd.format(**locals()), "Extract SV only regions for BubbleTree")
    return out_file

def _is_possible_loh(rec, params, somatic_info):
    """Check if the VCF record is a het in the normal with sufficient support.
    """
    normal_good, tumor_good = False, False
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
    return normal_good and tumor_good

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
            print rec
