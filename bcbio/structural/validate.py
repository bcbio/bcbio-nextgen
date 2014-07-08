"""Provide validation of structural variations against truth sets.

Tests overlaps of the combined ensemble structural variant BED against
a set of known regions.  Requires any overlap between ensemble set and
known regions, and removes regions from analysis that overlap with
exclusion regions.
"""

import toolz as tz
try:
    import pybedtools
except ImportError:
    pybedtools = None

def _stat_str(x, n):
    return "%.1f%% (%s / %s)" % (float(x) / float(n) * 100.0, x, n)

def _evaluate_one(caller, svtype, ensemble, truth, exclude):
    """Compare a ensemble results for a caller against a specific caller and SV type.
    """
    def is_caller_svtype(feat):
        for name in feat.name.split(","):
            if name.startswith((svtype, "cnv")) and (caller == "ensemble" or name.endswith(caller)):
                return True
        return False
    exfeats = pybedtools.BedTool(exclude)
    efeats = pybedtools.BedTool(ensemble).filter(is_caller_svtype).intersect(exfeats, v=True, f=0.25, r=True).saveas()
    tfeats = pybedtools.BedTool(truth).intersect(exfeats, v=True, f=0.25, r=True)
    etotal = len(set(efeats))
    ttotal = len(set(tfeats))
    match = len(set(efeats.intersect(tfeats)))
    return {"sensitivity": _stat_str(match, ttotal),
            "specificity": _stat_str(match, etotal)}

def _evaluate_multi(callers, truth_svtypes, ensemble, exclude):
    for svtype, truth in truth_svtypes.items():
        for caller in callers:
            evalout = _evaluate_one(caller, svtype, ensemble, truth, exclude)
            print caller, svtype, evalout

def evaluate(data, sv_calls):
    """Provide evaluations for multiple callers split by structural variant type.
    """
    truth_sets = tz.get_in(["config", "algorithm", "svvalidate"], data)
    print truth_sets
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
                    {"DEL": "NA12878.50X.ldgp.molpb_val.20140508.chr21.bed"},
                    "NA12878-ensemble.bed", "LCR.bed.gz")
