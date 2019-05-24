"""Variant quality control summaries

TODO:
 - Add plot of depth metrics to replace GATK based depth metrics calculation
"""
import os
import shutil

import six
import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.variation import genotype, germline

def run(_, data, out_dir):
    """Prepare variants QC analysis: bcftools stats and snpEff output.
    """
    out = []
    vcinfo = get_active_vcinfo(data)
    if vcinfo:
        if dd.get_phenotype(data) == "normal" and "germline" in vcinfo:
            out.append(_bcftools_stats(data, out_dir, "germline", germline=True))
        elif dd.get_phenotype(data) != "germline":
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

def _add_filename_details(full_f):
    """Add variant callers and germline information standard CWL filenames.

    This is an ugly way of working around not having metadata with calls.
    """
    out = {"vrn_file": full_f}
    f = os.path.basename(full_f)
    for vc in list(genotype.get_variantcallers().keys()) + ["ensemble"]:
        if f.find("-%s.vcf" % vc) > 0:
            out["variantcaller"] = vc
    if f.find("-germline-") >= 0:
        out["germline"] = full_f
    return out

def _get_variants(data):
    """Retrieve variants from CWL and standard inputs for organizing variants.
    """
    active_vs = []
    if "variants" in data:
        variants = data["variants"]
        # CWL based list of variants
        if isinstance(variants, dict) and "samples" in variants:
            variants = variants["samples"]
        for v in variants:
            # CWL -- a single variant file
            if isinstance(v, six.string_types) and os.path.exists(v):
                active_vs.append(_add_filename_details(v))
            elif (isinstance(v, (list, tuple)) and len(v) > 0 and
                  isinstance(v[0], six.string_types) and os.path.exists(v[0])):
                for subv in v:
                    active_vs.append(_add_filename_details(subv))
            elif isinstance(v, dict):
                if v.get("vrn_file"):
                    active_vs.append(v)
                elif v.get("population"):
                    vrnfile = v.get("population").get("vcf")
                    active_vs.append(_add_filename_details(vrnfile))
                elif v.get("vcf"):
                    active_vs.append(_add_filename_details(v.get("vcf")))
    return active_vs

def get_active_vcinfo(data, use_ensemble=True):
    """Use first caller if ensemble is not active
    """
    active_vs = _get_variants(data)
    if len(active_vs) > 0:
        e_active_vs = []
        if use_ensemble:
            e_active_vs = [v for v in active_vs if v.get("variantcaller") == "ensemble"]
        if len(e_active_vs) == 0:
            e_active_vs = [v for v in active_vs if v.get("variantcaller") != "ensemble"]
        if len(e_active_vs) > 0:
            return e_active_vs[0]

def extract_germline_vcinfo(data, out_dir):
    """Extract germline VCFs from existing tumor inputs.
    """
    supported_germline = set(["vardict", "octopus", "freebayes"])
    if dd.get_phenotype(data) in ["tumor"]:
        for v in _get_variants(data):
            if v.get("variantcaller") in supported_germline:
                if v.get("germline"):
                    return v
                else:
                    d = utils.deepish_copy(data)
                    d["vrn_file"] = v["vrn_file"]
                    gd = germline.extract(d, [d], out_dir)
                    v["germline"] = gd["vrn_file_plus"]["germline"]
                    return v
