"""Variant quality control summaries

TODO:
 - Add plot of depth metrics to replace GATK based depth metrics calculation
"""
import os
import shutil

import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.provenance import do

def run(bam_file, data, out_dir):
    """Prepare variants QC analysis: bcftools stats and snpEff output.
    """
    out = []
    for stat_fn in [_bcftools_stats, _snpeff_stats]:
        out_file = stat_fn(data, out_dir)
        if out_file:
            out.append(out_file)
    if out:
        return {"base": out[0], "secondary": out[1:]}

def _snpeff_stats(data, out_dir):
    vcinfo = get_active_vcinfo(data)
    if vcinfo and vcinfo.get("vrn_stats"):
        effects_csv = tz.get_in(["vrn_stats", "effects-stats-csv"], vcinfo)
        if effects_csv and utils.file_exists(effects_csv):
            out_dir = utils.safe_makedir(out_dir)
            out_file = os.path.join(out_dir, "%s.txt" % dd.get_sample_name(data))
            with file_transaction(data, out_file) as tx_out_file:
                shutil.copy(effects_csv, tx_out_file)
            return out_file

def get_active_vcinfo(data):
    """Use first caller if ensemble is not active
    """
    callers = dd.get_variantcaller(data)
    if not callers:
        return None
    if isinstance(callers, basestring):
        callers = [callers]
    active_vs = []
    if "variants" in data:
        for v in data["variants"]:
            if v.get("variantcaller") == "ensemble":
                return v
            if v.get("vrn_file"):
                active_vs.append(v)
        if len(active_vs) > 0:
            return active_vs[0]

def _bcftools_stats(data, out_dir):
    """Run bcftools stats.
    """
    vcinfo = get_active_vcinfo(data)
    if vcinfo:
        out_dir = utils.safe_makedir(out_dir)
        vcf_file = vcinfo["vrn_file"]
        if tz.get_in(("config", "algorithm", "jointcaller"), data):
            opts = ""
        else:
            opts = "-f PASS"
        name = dd.get_sample_name(data)
        out_file = os.path.join(out_dir, "%s_bcftools_stats.txt" % name)
        bcftools = config_utils.get_program("bcftools", data["config"])
        if not utils.file_exists(out_file):
            with file_transaction(data, out_file) as tx_out_file:
                orig_out_file = os.path.join(os.path.dirname(tx_out_file), "orig_%s" % os.path.basename(tx_out_file))
                cmd = ("{bcftools} stats -s {name} {opts} {vcf_file} > {orig_out_file}")
                do.run(cmd.format(**locals()), "bcftools stats %s" % dd.get_sample_name(data))
                with open(orig_out_file) as in_handle:
                    with open(tx_out_file, "w") as out_handle:
                        for line in in_handle:
                            if line.startswith("ID\t"):
                                parts = line.split("\t")
                                parts[-1] = "%s\n" % name
                                line = "\t".join(parts)
                            out_handle.write(line)
        return out_file
