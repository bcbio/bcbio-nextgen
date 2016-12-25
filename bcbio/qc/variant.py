"""Variant quality control summaries
"""
import os

import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.provenance import do

def run(bam_file, data, out_dir):
    """Prepare variants QC analysis: bcftools stats and snpEff output.
    """
    stats_file = bcftools_stats(data, out_dir)
    if stats_file:
        return {"base": stats_file, "secondary": []}

def _get_variant_callers(data):
    """Use first caller if ensemble is not active"""
    callers = dd.get_variantcaller(data)
    if not callers:
        return None
    if isinstance(callers, basestring):
        callers = [callers]
    active_callers = [c.get("variantcaller") for c in data.get("variants", [{}])]
    active_vcf = [c.get("vrn_file") for c in data.get("variants", [{}])]
    active_germline = [c.get("germline") for c in data.get("variants", [{}])]
    vcf = dict(zip(active_callers, active_vcf))
    germline = dict(zip(active_callers, active_germline))
    if "ensemble" in active_callers:
        vcf_fn = vcf["ensemble"]
    else:
        vcf_fn = vcf[callers[0]]
    if not vcf_fn:
        vcf_fn = germline[callers[0]]
    return vcf_fn

def bcftools_stats(data, out_dir):
    """Run bcftools stats.
    """
    vcf_file = _get_variant_callers(data)
    if vcf_file:
        if tz.get_in(("config", "algorithm", "jointcaller"), data):
            opts = ""
        else:
            opts = "-f PASS"
        name = dd.get_sample_name(data)
        stem = os.path.join(out_dir, os.path.basename(utils.splitext_plus(vcf_file)[0]))
        out_file = os.path.join(out_dir, "%s-%s-bcfstats.tsv" % (stem, name))
        bcftools = config_utils.get_program("bcftools", data["config"])
        if not utils.file_exists(out_file):
            with file_transaction(data, out_file) as tx_out_file:
                cmd = ("{bcftools} stats -s {name} {opts} {vcf_file} > {tx_out_file}")
                do.run(cmd.format(**locals()), "bcftools stats %s" % dd.get_sample_name(data))
        return out_file
