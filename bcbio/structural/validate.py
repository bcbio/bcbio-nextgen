"""Provide validation of structural variations against truth sets.

Tests overlaps of the combined ensemble structural variant BED against
a set of known regions.  Requires any overlap between ensemble set and
known regions, and removes regions from analysis that overlap with
exclusion regions.
"""
import csv
import os

import toolz as tz
import numpy as np
import pandas as pd
import pybedtools

from bcbio.log import logger
from bcbio import utils
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.structural import convert
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.variation import bedutils, vcfutils, ploidy, validateplot

mpl = utils.LazyImport("matplotlib")
plt = utils.LazyImport("matplotlib.pyplot")
sns = utils.LazyImport("seaborn")

# -- VCF based validation

def _evaluate_vcf(calls, truth_vcf, work_dir, data):
    out_file = os.path.join(work_dir, os.path.join("%s-sv-validate.csv" % dd.get_sample_name(data)))
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                writer = csv.writer(out_handle)
                writer.writerow(["sample", "caller", "vtype", "metric", "value"])
                for call in calls:
                    for stats in _validate_caller_vcf(call["vrn_file"], truth_vcf, dd.get_callable_regions(data),
                                                      call["variantcaller"], data):

                        writer.writerow(stats)
    return out_file

def _validate_caller_vcf(call_vcf, truth_vcf, callable_regions, svcaller, data):
    """Validate a caller VCF against truth within callable regions, returning stratified stats
    """
    stats = _calculate_comparison_stats(truth_vcf)
    callable_regions = bedutils.clean_file(callable_regions, data)
    callable_bed = pybedtools.BedTool(callable_regions).merge(d=stats["merge_size"]).saveas().fn

    match_calls = set([])
    truth_stats = {"tp": [], "fn": [], "fp": []}
    calls_by_region = {}
    call_vcf = slim_vcf(call_vcf, data)
    for call in _callable_intersect(call_vcf, callable_bed, data):
        key = tuple(call[:5] + call[7:8])
        calls_by_region[tuple(call[-3:])] = key

    truth = None
    regions = []
    for parts in _callable_intersect(truth_vcf, callable_bed, data):
        cur_region = tuple(parts[-3:])
        cur_truth = tuple(parts[:5] + parts[7:8])
        if truth is None:
            truth = cur_truth
        if cur_truth == truth:
            regions.append(cur_region)
        else:
            match_calls, truth_stats = _check_call(truth, regions, calls_by_region, match_calls, truth_stats)
            truth = cur_truth
            regions = [cur_region]
    with utils.open_gzipsafe(call_vcf) as in_handle:
        for call in (l.split("\t") for l in in_handle if not l.startswith("#")):
            start, end = _get_start_end(call)
            if end:
                key = tuple(call[:5] + call[7:8])
                if key not in match_calls:
                    call_info = _summarize_call(key)
                    if _event_passes(call_info, stats):
                        truth_stats["fp"].append(call_info)
    return _to_csv(truth_stats, stats, dd.get_sample_name(data), svcaller)

def _to_csv(truth_stats, stats, sample, svcaller):
    out = []
    for metric, vals in truth_stats.items():
        for svtype in sorted(list(stats["svtypes"])):
            count = len([x for x in vals if x["svtype"] == svtype])
            out.append([sample, svcaller, svtype, metric, count])
            for start, end in stats["ranges"]:
                count = len([x for x in vals if (x["svtype"] == svtype
                                                 and x["size"] >= start and x["size"] < end)])
                out.append([sample, svcaller, "%s_%s-%s" % (svtype, start, end), metric, count])
    return out

def _event_passes(call_info, stats):
    return (call_info["svtype"] in stats["svtypes"] and call_info["size"] > stats["min_size"]
            and call_info["size"] < stats["max_size"])

def _calculate_comparison_stats(truth_vcf):
    """Identify calls to validate from the input truth VCF.
    """
    sizes = []
    svtypes = set([])
    with utils.open_gzipsafe(truth_vcf) as in_handle:
        for call in (l.rstrip().split("\t") for l in in_handle if not l.startswith("#")):
            stats = _summarize_call(tuple(call[:5] + call[7:8]))
            sizes.append(stats["size"])
            svtypes.add(stats["svtype"])
    pct10 = int(np.percentile(sizes, 10))
    pct25 = int(np.percentile(sizes, 25))
    pct50 = int(np.percentile(sizes, 50))
    pct75 = int(np.percentile(sizes, 75))
    ranges_detailed = [(int(min(sizes)), pct10), (pct10, pct25), (pct25, pct50),
                       (pct50, pct75), (pct75, max(sizes))]
    ranges_split = [(int(min(sizes)), pct50), (pct50, max(sizes))]
    return {"min_size": int(min(sizes) * 0.95), "max_size": int(max(sizes) + 1.05),
            "svtypes": svtypes, "merge_size": int(np.percentile(sizes, 10)),
            "ranges": []}

def _get_start_end(parts, index=7):
    """Retrieve start and end for a VCF record, skips BNDs without END coords
    """
    start = parts[1]
    end = [x.split("=")[-1] for x in parts[index].split(";") if x.startswith("END=")]
    if end:
        end = end[0]
        return start, end
    return None, None

def _callable_intersect(in_file, callable_bed, data):
    """Return list of original VCF SVs intersected by callable regions.

    Does not try to handle BNDs. We should resolve these and return where possible.
    """
    with tx_tmpdir(data) as tmpdir:
        in_bed = os.path.join(tmpdir, "%s-convert.bed" % utils.splitext_plus(os.path.basename(in_file))[0])
        with utils.open_gzipsafe(in_file) as in_handle:
            with open(in_bed, "w") as out_handle:
                for parts in (l.split("\t") for l in in_handle if not l.startswith("#")):
                    start, end = _get_start_end(parts)
                    if end:
                        out_handle.write("\t".join([parts[0], start, end] + parts) + "\n")
        out_file = os.path.join(tmpdir, "%s-subset.tsv" % utils.splitext_plus(os.path.basename(in_file))[0])
        cmd = "bedtools intersect -a {in_bed} -b {callable_bed} -wa -wb > {out_file}"
        do.run(cmd.format(**locals()), "Intersect VCF by callable")
        with open(out_file) as in_handle:
            for line in in_handle:
                yield line.rstrip().split("\t")[3:]

def _summarize_call(parts):
    """Provide summary metrics on size and svtype for a SV call.
    """
    svtype = [x.split("=")[-1] for x in parts[-1].split(";") if x.startswith("SVTYPE=")]
    svtype = svtype[0] if svtype else ""
    start, end = _get_start_end(parts, index=-1)
    return {"svtype": svtype, "size": int(end) - int(start)}

def _is_compatible(call, truth):
    """Determine if a call and truth set are equivalent matches.

    Entry point to restrict matches by call type or size.
    """
    call = _summarize_call(call)
    truth = _summarize_call(truth)
    return call["size"] < (truth["size"] * 2)

def _check_call(truth, regions, calls_by_region, match_calls, truth_stats):
    """Validate if we have a call that matches this given truth set within callable regions.
    """
    cur_matches = []
    for region in regions:
        cur_match_call = calls_by_region.get(region)
        if cur_match_call:
            if cur_match_call not in match_calls and _is_compatible(cur_match_call, truth):
                cur_matches.append(cur_match_call)
            match_calls.add(cur_match_call)
    if cur_matches:
        truth_stats["tp"].append(_summarize_call(truth))
    else:
        truth_stats["fn"].append(_summarize_call(truth))
    return match_calls, truth_stats

def slim_vcf(in_file, data):
    """Remove larger annotations which slow down VCF processing
    """
    to_remove = ["ANN", "LOF"]
    to_remove_str = tuple(["##INFO=<ID=%s" % x for x in to_remove])
    in_file = vcfutils.bgzip_and_index(in_file, data, remove_orig=False)
    out_file = "%s-slim.vcf.gz" % utils.splitext_plus(in_file)[0]
    if not utils.file_uptodate(out_file, in_file):
        cur_remove = []
        with utils.open_gzipsafe(in_file) as in_handle:
            for line in in_handle:
                if not line.startswith("#"):
                    break
                elif line.startswith(to_remove_str):
                    cur_id = line.split("ID=")[-1].split(",")[0]
                    cur_remove.append("INFO/%s" % cur_id)
        with file_transaction(data, out_file) as tx_out_file:
            if cur_remove:
                cur_remove = ",".join(cur_remove)
                cmd = ("bcftools view -f 'PASS,.' {in_file} | "
                       "bcftools annotate -x {cur_remove} -O z -o {tx_out_file}")
            else:
                cmd = ("bcftools view -f 'PASS,.' {in_file} -O z -o {tx_out_file}")
            do.run(cmd.format(**locals()), "Create slim VCF")
    return out_file

# -- BED based evaluation

EVENT_SIZES = [(100, 450), (450, 2000), (2000, 4000), (4000, 20000), (20000, 60000),
               (60000, int(1e6))]

def _stat_str(x, n):
    if n > 0:
        val = float(x) / float(n) * 100.0
        return {"label": "%.1f%% (%s / %s)" % (val, x, n), "val": val}
    else:
        return {"label": "", "val": 0}

def cnv_to_event(name, data):
    """Convert a CNV to an event name.
    """
    cur_ploidy = ploidy.get_ploidy([data])
    if name.startswith("cnv"):
        num = max([int(x) for x in name.split("_")[0].replace("cnv", "").split(";")])
        if num < cur_ploidy:
            return "DEL"
        elif num > cur_ploidy:
            return "DUP"
        else:
            return name
    else:
        return name

def _evaluate_one(caller, svtype, size_range, ensemble, truth, data):
    """Compare a ensemble results for a caller against a specific caller and SV type.
    """
    def cnv_matches(name):
        return cnv_to_event(name, data) == svtype
    def is_breakend(name):
        return name.startswith("BND")
    def in_size_range(max_buffer=0):
        def _work(feat):
            minf, maxf = size_range
            buffer = min(max_buffer, int(((maxf + minf) / 2.0) / 10.0))
            size = feat.end - feat.start
            return size >= max([0, minf - buffer]) and size < maxf + buffer
        return _work
    def is_caller_svtype(feat):
        for name in feat.name.split(","):
            if ((name.startswith(svtype) or cnv_matches(name) or is_breakend(name))
                  and (caller == "sv-ensemble" or name.endswith(caller))):
                return True
        return False
    minf, maxf = size_range
    efeats = pybedtools.BedTool(ensemble).filter(in_size_range(0)).filter(is_caller_svtype).saveas().sort().merge()
    tfeats = pybedtools.BedTool(truth).filter(in_size_range(0)).sort().merge().saveas()
    etotal = efeats.count()
    ttotal = tfeats.count()
    match = efeats.intersect(tfeats, u=True).sort().merge().saveas().count()
    return {"sensitivity": _stat_str(match, ttotal),
            "precision": _stat_str(match, etotal)}

def _evaluate_multi(calls, truth_svtypes, work_dir, data):
    base = os.path.join(work_dir, "%s-sv-validate" % (dd.get_sample_name(data)))
    out_file = base + ".csv"
    df_file = base + "-df.csv"
    if any((not utils.file_uptodate(out_file, x["vrn_file"])
            or not utils.file_uptodate(df_file, x["vrn_file"])) for x in calls):
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                with open(df_file, "w") as df_out_handle:
                    writer = csv.writer(out_handle)
                    dfwriter = csv.writer(df_out_handle)
                    writer.writerow(["svtype", "size", "caller", "sensitivity", "precision"])
                    dfwriter.writerow(["svtype", "size", "caller", "metric", "value", "label"])
                    for svtype, truth in truth_svtypes.items():
                        for size in EVENT_SIZES:
                            str_size = "%s-%s" % size
                            for call in calls:
                                call_bed = convert.to_bed(call, dd.get_sample_name(data), work_dir, calls, data)
                                if utils.file_exists(call_bed):
                                    evalout = _evaluate_one(call["variantcaller"], svtype, size, call_bed,
                                                            truth, data)
                                    writer.writerow([svtype, str_size, call["variantcaller"],
                                                     evalout["sensitivity"]["label"], evalout["precision"]["label"]])
                                    for metric in ["sensitivity", "precision"]:
                                        dfwriter.writerow([svtype, str_size, call["variantcaller"], metric,
                                                           evalout[metric]["val"], evalout[metric]["label"]])
    return out_file, df_file

def _plot_evaluation(df_csv):
    if mpl is None or plt is None or sns is None:
        not_found = ", ".join([x for x in ['mpl', 'plt', 'sns'] if eval(x) is None])
        logger.info("No validation plot. Missing imports: %s" % not_found)
        return None
    mpl.use('Agg', force=True)
    df = pd.read_csv(df_csv).fillna("0%")
    out = {}
    for event in df["svtype"].unique():
        out[event] = _plot_evaluation_event(df_csv, event)
    return out

def _plot_evaluation_event(df_csv, svtype):
    """Provide plot of evaluation metrics for an SV event, stratified by event size.
    """
    titles = {"INV": "Inversions", "DEL": "Deletions", "DUP": "Duplications",
              "INS": "Insertions"}
    out_file = "%s-%s.png" % (os.path.splitext(df_csv)[0], svtype)
    sns.set(style='white')
    if not utils.file_uptodate(out_file, df_csv):
        metrics = ["sensitivity", "precision"]
        df = pd.read_csv(df_csv).fillna("0%")
        df = df[(df["svtype"] == svtype)]
        event_sizes = _find_events_to_include(df, EVENT_SIZES)
        fig, axs = plt.subplots(len(event_sizes), len(metrics), tight_layout=True)
        if len(event_sizes) == 1:
            axs = [axs]
        callers = sorted(df["caller"].unique())
        if "sv-ensemble" in callers:
            callers.remove("sv-ensemble")
            callers.append("sv-ensemble")
        for i, size in enumerate(event_sizes):
            size_label = "%s to %sbp" % size
            size = "%s-%s" % size
            for j, metric in enumerate(metrics):
                ax = axs[i][j]
                ax.get_xaxis().set_ticks([])
                ax.spines['bottom'].set_visible(False)
                ax.spines['left'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.set_xlim(0, 125.0)
                if i == 0:
                    ax.set_title(metric, size=12, y=1.2)
                vals, labels = _get_plot_val_labels(df, size, metric, callers)
                ax.barh(np.arange(len(vals)), vals)
                if j == 0:
                    ax.tick_params(axis='y', which='major', labelsize=8)
                    ax.locator_params(nbins=len(callers) + 2, axis="y", tight=True)
                    ax.set_yticklabels(callers, va="bottom")
                    ax.text(100, len(callers), size_label, fontsize=10)
                else:
                    ax.get_yaxis().set_ticks([])
                for ai, (val, label) in enumerate(zip(vals, labels)):
                    ax.annotate(label, (val + 0.75, ai + 0.35), va='center', size=7)
        if svtype in titles:
            fig.text(0.025, 0.95, titles[svtype], size=14)
        fig.set_size_inches(7, len(event_sizes) + 1)
        fig.savefig(out_file)
    return out_file

def _find_events_to_include(df, event_sizes):
    out = []
    for size in event_sizes:
        str_size = "%s-%s" % size
        curdf = df[(df["size"] == str_size) & (df["metric"] == "sensitivity")]
        for val in list(curdf["label"]):
            if val != "0%":
                out.append(size)
                break
    return out

def _get_plot_val_labels(df, size, metric, callers):
    curdf = df[(df["size"] == size) & (df["metric"] == metric)]
    vals, labels = [], []
    for caller in callers:
        row = curdf[curdf["caller"] == caller]
        val = list(row["value"])[0]
        if val == 0:
            val = 0.1
        vals.append(val)
        labels.append(list(row["label"])[0])
    return vals, labels

# -- general functionality

def evaluate(data):
    """Provide evaluations for multiple callers split by structural variant type.
    """
    work_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                               dd.get_sample_name(data), "validate"))
    truth_sets = tz.get_in(["config", "algorithm", "svvalidate"], data)
    if truth_sets and data.get("sv"):
        if isinstance(truth_sets, dict):
            val_summary, df_csv = _evaluate_multi(data["sv"], truth_sets, work_dir, data)
            summary_plots = _plot_evaluation(df_csv)
            data["sv-validate"] = {"csv": val_summary, "plot": summary_plots, "df": df_csv}
        else:
            assert isinstance(truth_sets, basestring) and utils.file_exists(truth_sets), truth_sets
            val_summary = _evaluate_vcf(data["sv"], truth_sets, work_dir, data)
            title = "%s structural variants" % dd.get_sample_name(data)
            summary_plots = validateplot.classifyplot_from_valfile(val_summary, outtype="png", title=title)
            data["sv-validate"] = {"csv": val_summary, "plot": summary_plots[0]}
    return data

if __name__ == "__main__":
    #_, df_csv = _evaluate_multi(["lumpy", "delly", "wham", "sv-ensemble"],
    #                            {"DEL": "synthetic_challenge_set3_tumor_20pctmasked_truth_sv_DEL.bed"},
    #                            "syn3-tumor-ensemble-filter.bed", "sv_exclude.bed")
    #_, df_csv = _evaluate_multi(["lumpy", "delly", "cn_mops", "sv-ensemble"],
    #                            {"DEL": "NA12878.50X.ldgp.molpb_val.20140508.bed"},
    #                            "NA12878-ensemble.bed", "LCR.bed.gz")
    import sys
    _plot_evaluation(sys.argv[1])
