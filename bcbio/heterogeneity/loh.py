"""Summarize loss of heterozygosity (LOH) from heterogeneity callers.

Provides high level summaries of LOH calls in regions of interest.
"""
import collections
import os

import yaml

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd

_COORDS = {"hg38": {"HLA": ("chr6", 28510120, 33480577),
                    "B2M": ("chr15", 44711487, 44718877),
                    "PTEN": ("chr10", 87863113, 87971930),
                    "CDKN2A": ("chr9", 21967753, 21995301)},
           "hg19": {"HLA": ("chr6", 29640000, 33120000),
                    "B2M": ("chr15", 45003675, 45011075),
                    "PTEN": ("chr10", 89622870, 89731687),
                    "CDKN2A": ("chr9", 21967751, 21995300)},
           "GRCh37": {"HLA": ("6", 29640000, 33120000),
                      "B2M": ("15", 45003675, 45011075),
                      "PTEN": ("10", 89622870, 89731687),
                      "CDKN2A": ("9", 21967751, 21995300)}}

def get_coords(data):
    return _COORDS.get(dd.get_genome_build(data) or {})

def summary_status(call, data):
    """Retrieve LOH status in regions of interest, along with heterogeneity metrics.

    Provides output with overall purity and ploidy, along with region
    specific LOH.
    """
    coords = get_coords(data)
    out_file = None
    if call.get("vrn_file") and os.path.exists(call.get("vrn_file")):
        out_file = os.path.join(os.path.dirname(call["vrn_file"]),
                                "%s-%s-lohsummary.yaml" % (dd.get_sample_name(data), call["variantcaller"]))
        if not utils.file_uptodate(out_file, call["vrn_file"]):
            out = {}
            if call["variantcaller"] == "titancna":
                out.update(_titancna_summary(call, coords))
                pass
            elif call["variantcaller"] == "purecn":
                out.update(_purecn_summary(call, coords))
            if out:
                out["description"] = dd.get_sample_name(data)
                out["variantcaller"] = call["variantcaller"]
                with file_transaction(data, out_file) as tx_out_file:
                    with open(tx_out_file, "w") as out_handle:
                        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out_file if out_file and os.path.exists(out_file) else None

def _titancna_summary(call, coords):
    """Summarize purity, ploidy and LOH for TitanCNA.
    """
    out = {}
    loh_calls = {k: collections.defaultdict(int) for k in coords.keys()}
    with open(call["subclones"]) as in_handle:
        header = in_handle.readline().strip().split()
        for line in in_handle:
            val = dict(zip(header, line.strip().split()))
            start = int(val["Start_Position.bp."])
            end = int(val["End_Position.bp."])
            for region, cur_coords in coords.items():
                if val["Chromosome"] == cur_coords[0] and are_overlapping((start, end), cur_coords[1:]):
                    if int(val["MinorCN"]) == 0:
                        loh_calls[region]["loh"] += 1
                    else:
                        loh_calls[region]["std"] += 1
    out["LOH"] = {r: _get_loh_from_calls(c) for r, c in loh_calls.items()}
    with open(call["hetsummary"]) as in_handle:
        vals = dict(zip(in_handle.readline().strip().split("\t"), in_handle.readline().strip().split("\t")))
    out["purity"] = vals["purity"]
    out["ploidy"] = vals["ploidy"]
    return out

def _purecn_summary(call, coords):
    """Summarize purity, ploidy and LOH for PureCN.
    """
    out = {}
    loh_calls = {k: collections.defaultdict(int) for k in coords.keys()}
    with open(call["loh"]) as in_handle:
        in_handle.readline()  # header
        for line in in_handle:
            _, chrom, start, end, _, cn, minor_cn = line.split(",")[:7]
            start = int(start)
            end = int(end)
            for region, cur_coords in coords.items():
                if chrom == cur_coords[0] and are_overlapping((start, end), cur_coords[1:]):
                    if int(minor_cn) == 0:
                        loh_calls[region]["loh"] += 1
                    else:
                        loh_calls[region]["std"] += 1
    out["LOH"] = {r: _get_loh_from_calls(c) for r, c in loh_calls.items()}
    with open(call["hetsummary"]) as in_handle:
        vals = dict(zip(in_handle.readline().strip().replace('"', '').split(","),
                        in_handle.readline().strip().split(",")))
    out["purity"] = vals["Purity"]
    out["ploidy"] = vals["Ploidy"]
    return out

def _get_loh_from_calls(calls):
    if calls["loh"]:
        return "mixed LOH" if calls["std"] else "LOH"
    else:
        return "no LOH"

def are_overlapping(r, s):
    """Test if two coordinates overlap.
    https://stackoverflow.com/a/27182551
    """
    return r[1] >= s[0] and s[1] >= r[0]
