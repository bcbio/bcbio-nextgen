"""Estimate tumor purity and frequency using copy number and allele freqencies with BubbleTree.

http://www.bioconductor.org/packages/release/bioc/html/BubbleTree.html
http://www.bioconductor.org/packages/release/bioc/vignettes/BubbleTree/inst/doc/BubbleTree-vignette.html
"""
import os

from pysam import VariantFile

from bcbio import utils
from bcbio.distributed.transaction import file_transaction

def run(vrn_info, cnv_info, somatic_info):
    """Run BubbleTree given variant calls, CNVs and somatic
    """
    work_dir = _cur_workdir(somatic_info.tumor_data)
    input_vcf = _prep_vrn_file(vrn_info["vrn_file"], work_dir, somatic_info)
    print input_vcf

def _prep_vrn_file(in_file, work_dir, somatic_info):
    """Select heterozygous variants in the normal sample with sufficient depth.
    """
    data = somatic_info.tumor_data
    params = {"min_freq": 0.4,
              "max_freq": 0.6,
              "min_depth": 15}
    out_file = os.path.join(work_dir, "%s-prep.vcf" % (utils.splitext_plus(os.path.basename(in_file))))
    if not utils.file_uptodata(out_file, in_file):
        with file_transaction(data, out_file) as tx_out_file:
            bcf_in = VariantFile(in_file)
            bcf_out = VariantFile(tx_out_file, "w", header=bcf_in.header)
            for rec in bcf_in:
                if _is_possible_loh(rec, params, somatic_info):
                    bcf_out.write(rec)
    return out_file

def _is_possible_loh(rec, params, somatic_info):
    """Check if the VCF record is a het in the normal with sufficient support.
    """
    normal_good, tumor_good = False, False
    for sample in rec.samples:
        if sample.name == somatic_info.normal_name:
            pass
        elif sample.name == somatic_info.tumor_name:
            pass
    return normal_good and tumor_good

def _cur_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "heterogeneity",
                                           dd.get_sample_name(data), "bubbletree"))