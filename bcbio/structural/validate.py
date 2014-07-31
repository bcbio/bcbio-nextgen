"""Provide validation of structural variations against truth sets.

Tests overlaps of the combined ensemble structural variant BED against
a set of known regions.  Requires any overlap between ensemble set and
known regions, and removes regions from analysis that overlap with
exclusion regions.
"""
import collections
import csv

import toolz as tz
try:
    import pybedtools
except ImportError:
    pybedtools = None

from bcbio import utils

def _stat_str(x, n):
    if n > 0:
        return "%.1f%% (%s / %s)" % (float(x) / float(n) * 100.0, x, n)
    else:
        return "0"

def _evaluate_one(caller, svtype, size_range, ensemble, truth, exclude):
    """Compare a ensemble results for a caller against a specific caller and SV type.
    """
    def cnv_matches(name):
        """Check for CNV matches -- XXX hardcoded for diploid comparisons.
        """
        if name.startswith("cnv"):
            num = int(name.split("_")[0].replace("cnv", ""))
            if svtype == "DEL" and num < 2:
                return True
            elif svtype == "DUP" and num > 2:
                return True
            else:
                return False
        else:
            return False
    def in_size_range(feat):
        minf, maxf = size_range
        size = feat.end - feat.start
        return size >= minf and size < maxf
    def is_caller_svtype(feat):
        for name in feat.name.split(","):
            if (name.startswith(svtype) or cnv_matches(name)) and (caller == "ensemble" or name.endswith(caller)):
                return True
        return False
    exfeats = pybedtools.BedTool(exclude)
    efeats = pybedtools.BedTool(ensemble).filter(in_size_range).filter(is_caller_svtype).saveas()\
                       .intersect(exfeats, v=True, f=0.50, r=True).sort().merge().saveas()
    tfeats = pybedtools.BedTool(truth).filter(in_size_range)\
                                      .intersect(exfeats, v=True, f=0.50, r=True).sort().merge().saveas()
    etotal = efeats.count()
    ttotal = tfeats.count()
    match = efeats.intersect(tfeats).sort().merge().saveas().count()
    return {"sensitivity": _stat_str(match, ttotal),
            "precision": _stat_str(match, etotal)}

def _evaluate_multi(callers, truth_svtypes, ensemble, exclude):
    #_validated_stats(truth_svtypes, ensemble, exclude)
    out_file = "%s-validate.csv" % utils.splitext_plus(ensemble)[0]
    sizes = [(1, 250), (250, 1000), (1000, 5000), (5000, 25000), (25000, int(1e6))]
    if not utils.file_uptodate(out_file, ensemble):
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle)
            writer.writerow(["svtype", "size", "caller", "sensitivity", "precision"])
            for svtype, truth in truth_svtypes.items():
                for size in sizes:
                    for caller in callers:
                        evalout = _evaluate_one(caller, svtype, size, ensemble, truth, exclude)
                        writer.writerow([svtype, "%s-%s" % size, caller,
                                         evalout["sensitivity"], evalout["precision"]])
    return out_file

def _validated_stats(truth_svtypes, ensemble, exclude):
    """Provide statistics on which callers identify known structural variant calls.
    """
    for svtype, truth in truth_svtypes.items():
        exfeats = pybedtools.BedTool(exclude)
        tfeats = pybedtools.BedTool(truth).intersect(exfeats, v=True, f=0.50, r=True).sort().merge().saveas()
        overlaps = collections.defaultdict(set)
        for feat in tfeats.intersect(ensemble, wao=True):
            size = feat.end - feat.start
            callers, overlap = feat.fields[-2:]
            if int(overlap) > 0:
                callers = set([x.split("_", 1)[-1] for x in callers.split(",")])
            else:
                callers = set([])
            key = (size, feat.chrom, feat.start, feat.end)
            overlaps[key] = overlaps[key].union(callers)
        for key in sorted(overlaps.keys()):
            print key[0], sorted(list(overlaps[key]))

def evaluate(data, sv_calls):
    """Provide evaluations for multiple callers split by structural variant type.
    """
    truth_sets = tz.get_in(["config", "algorithm", "svvalidate"], data)
    ensemble_callsets = [(i, x) for (i, x) in enumerate(sv_calls) if x["variantcaller"] == "ensemble"]
    if truth_sets and len(ensemble_callsets) > 0:
        callers = [x["variantcaller"] for x in sv_calls]
        ensemble_bed = ensemble_callsets[0][-1]["vrn_file"]
        exclude_files = [f for f in [x.get("exclude_file") for x in sv_calls] if f]
        exclude_file = exclude_files[0] if len(exclude_files) > 0 else None
        val_summary = _evaluate_multi(callers, truth_sets, ensemble_bed, exclude_file)
        ensemble_i, ensemble = ensemble_callsets[0]
        ensemble["validate"] = val_summary
        sv_calls[ensemble_i] = ensemble
    return sv_calls

if __name__ == "__main__":
    _evaluate_multi(["lumpy", "delly", "cn_mops", "ensemble"],
                    {"DEL": "NA12878.50X.ldgp.molpb_val.20140508.bed"},
                    "NA12878-ensemble.bed", "LCR.bed.gz")
