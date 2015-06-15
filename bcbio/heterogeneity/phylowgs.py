"""Reconstruct subclones and phylogenetic history using PhyloWGS.

PhyloWGS uses phylogenetic histories inferred through Battenberg CNV calls
along with variant frequencies.

https://github.com/morrislab/phylowgs
http://genomebiology.com/2015/16/1/35
"""
import os
import sys

from pysam import VariantFile

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do

def run(vrn_info, cnvs_by_name, somatic_info):
    """Run PhyloWGS given variant calls, CNVs and tumor/normal information.
    """
    work_dir = _cur_workdir(somatic_info.tumor_data)
    assert "battenberg" in cnvs_by_name, "PhyloWGS requires Battenberg CNV calls"
    ssm_file, cnv_file = _prep_inputs(vrn_info, cnvs_by_name["battenberg"], somatic_info, work_dir)
    print ssm_file, cnv_file

def _prep_inputs(vrn_info, cnv_info, somatic_info, work_dir):
    """Prepare inputs for running PhyloWGS from variant and CNV calls.
    """
    exe = os.path.join(os.path.dirname(sys.executable), "create_phylowgs_inputs.py")
    assert os.path.exists(exe), "Could not find input prep script for PhyloWGS runs."
    ssm_file = os.path.join(work_dir, "ssm_data.txt")
    cnv_file = os.path.join(work_dir, "cnv_data.txt")
    if not utils.file_exists(ssm_file) or not utils.file_exists(cnv_file):
        with file_transaction(somatic_info.tumor_data, ssm_file, cnv_file) as (tx_ssm_file, tx_cnv_file):
            variant_type, input_vcf_file = _prep_vrn_file(vrn_info["vrn_file"], vrn_info["variantcaller"],
                                                          work_dir, somatic_info)
            cmd = [sys.executable, exe,
                   "--battenberg", cnv_info["subclones"], "--contamination", _read_contam(cnv_info["contamination"]),
                   "--output-cnvs", tx_cnv_file, "--output-variants", tx_ssm_file,
                   "--variant-type", variant_type, input_vcf_file]
            do.run(cmd, "Prepare PhyloWGS inputs.")
    return ssm_file, cnv_file

def _prep_vrn_file(in_file, vcaller, work_dir, somatic_info):
    """Create a variant file to feed into the PhyloWGS prep script, limiting records.

    Handles MuTect and VarDict as inputs to PhyloWGS.
    """
    if vcaller.startswith("vardict"):
        variant_type = "vardict"
    elif vcaller == "mutect":
        variant_type = "mutect-smchet"
    else:
        raise ValueError("Unexpected variant caller for PhyloWGS prep: %s" % vcaller)
    out_file = os.path.join(work_dir, "%s-%s-prep.csv" % (utils.splitext_plus(os.path.basename(in_file))[0],
                                                          vcaller))
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(somatic_info.tumor_data, out_file) as tx_out_file:
            bcf_in = VariantFile(in_file)
            bcf_out = VariantFile(tx_out_file, "w", header=bcf_in.header)
            for rec in bcf_in:
                if _is_snp(rec) and "PASS" in rec.filter.keys():
                    bcf_out.write(rec)
    return variant_type, out_file

def _is_snp(rec):
    return max([len(x) for x in rec.alleles]) == 1

def _read_contam(in_file):
    with open(in_file) as in_handle:
        return in_handle.readline().strip()

def _cur_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "heterogeneity",
                                           dd.get_sample_name(data), "phylowgs"))
