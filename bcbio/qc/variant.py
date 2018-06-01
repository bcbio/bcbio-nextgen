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

def run(_, data, out_dir):
    """Prepare variants QC analysis: bcftools stats and snpEff output.
    """
    out = []
    vcinfo = get_active_vcinfo(data)
    if vcinfo:
        if dd.get_phenotype(data) != "germline":
            out.append(_bcftools_stats(data, out_dir))
            if "germline" in vcinfo:
                out.append(_bcftools_stats(data, out_dir, "germline", germline=True))
        else:
            out.append(_bcftools_stats(data, out_dir, germline=True))

    out.append(_snpeff_stats(data, out_dir))

    out = [item for item in out if item]
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

def _bcftools_stats(data, out_dir, vcf_file_key=None, germline=False):
    """Run bcftools stats.
    """
    vcinfo = get_active_vcinfo(data)
    if vcinfo:
        out_dir = utils.safe_makedir(out_dir)
        vcf_file = vcinfo[vcf_file_key or "vrn_file"]
        if dd.get_jointcaller(data) or "gvcf" in dd.get_tools_on(data):
            opts = ""
        else:
            opts = "-f PASS,."
        name = dd.get_sample_name(data)
        out_file = os.path.join(out_dir, "%s_bcftools_stats%s.txt" % (name, ("_germline" if germline else "")))
        bcftools = config_utils.get_program("bcftools", data["config"])
        if not utils.file_exists(out_file):
            with file_transaction(data, out_file) as tx_out_file:
                orig_out_file = os.path.join(os.path.dirname(tx_out_file), "orig_%s" % os.path.basename(tx_out_file))
                cmd = ("{bcftools} stats -s {name} {opts} {vcf_file} > {orig_out_file}")
                do.run(cmd.format(**locals()), "bcftools stats %s" % name)
                with open(orig_out_file) as in_handle:
                    with open(tx_out_file, "w") as out_handle:
                        for line in in_handle:
                            if line.startswith("ID\t"):
                                parts = line.split("\t")
                                parts[-1] = "%s\n" % name
                                line = "\t".join(parts)
                            out_handle.write(line)
        return out_file

def get_active_vcinfo(data, use_ensemble=True):
    """Use first caller if ensemble is not active

    Handles both CWL and standard inputs for organizing variants.
    """
    active_vs = []
    if "variants" in data:
        variants = data["variants"]
        # CWL based list of variants
        if isinstance(variants, dict) and "samples" in variants:
            variants = variants["samples"]
        for v in variants:
            # CWL -- a single variant file
            if isinstance(v, basestring) and os.path.exists(v):
                active_vs.append({"vrn_file": v})
            elif (isinstance(v, (list, tuple)) and len(v) > 0 and
                  isinstance(v[0], basestring) and os.path.exists(v[0])):
                active_vs.append({"vrn_file": v[0]})
            elif isinstance(v, dict) and v.get("variantcaller") == "ensemble":
                if use_ensemble:
                    return v
            elif isinstance(v, dict) and v.get("vrn_file"):
                active_vs.append(v)
        if len(active_vs) > 0:
            return active_vs[0]
