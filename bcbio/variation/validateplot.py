"""Plot validation results from variant calling comparisons.

Handles data normalization and plotting, emphasizing comparisons on methodology
differences.
"""
import collections
import os

import numpy as np
import pandas as pd

from bcbio.log import logger
from bcbio import utils
from bcbio.variation import bamprep

mpl = utils.LazyImport("matplotlib")
plt = utils.LazyImport("matplotlib.pyplot")
mpl_ticker = utils.LazyImport("matplotlib.ticker")
sns = utils.LazyImport("seaborn")

def classifyplot_from_plotfiles(plot_files, out_csv, outtype="png", title=None, size=None):
    """Create a plot from individual summary csv files with classification metrics.
    """
    dfs = [pd.read_csv(x) for x in plot_files]
    samples = []
    for df in dfs:
        for sample in df["sample"].unique():
            if sample not in samples:
                samples.append(sample)
    df = pd.concat(dfs)
    df.to_csv(out_csv, index=False)
    return classifyplot_from_valfile(out_csv, outtype, title, size, samples)

def classifyplot_from_valfile(val_file, outtype="png", title=None, size=None,
                              samples=None, callers=None):
    """Create a plot from a summarized validation file.

    Does new-style plotting of summarized metrics of
    false negative rate and false discovery rate.
    https://en.wikipedia.org/wiki/Sensitivity_and_specificity
    """
    mpl.use('Agg', force=True)
    df = pd.read_csv(val_file)
    grouped = df.groupby(["sample", "caller", "vtype"])
    df = grouped.apply(_calculate_fnr_fdr)
    df = df.reset_index()
    out_file = "%s.%s" % (os.path.splitext(val_file)[0], outtype)
    _do_classifyplot(df, out_file, title, size, samples, callers)
    return [out_file]

def _calculate_fnr_fdr(group):
    """Calculate the false negative rate (1 - sensitivity) and false discovery rate (1 - precision).
    """
    data = {k: d["value"] for k, d in group.set_index("metric").T.to_dict().items()}
    return pd.DataFrame([{"fnr": data["fn"] / float(data["tp"] + data["fn"]) * 100.0 if data["tp"] > 0 else 0.0,
                          "fdr": data["fp"] / float(data["tp"] + data["fp"]) * 100.0 if data["tp"] > 0 else 0.0,
                          "tpr": "TP: %s FN: %s" % (data["tp"], data["fn"]),
                          "spc": "FP: %s" % (data["fp"])}])

def _do_classifyplot(df, out_file, title=None, size=None, samples=None, callers=None):
    """Plot using classification-based plot using seaborn.
    """
    metric_labels = {"fdr": "False discovery rate",
                     "fnr": "False negative rate"}
    metrics = [("fnr", "tpr"), ("fdr", "spc")]
    colors = ["light grey", "greyish"]
    data_dict = df.set_index(["sample", "caller", "vtype"]).T.to_dict()
    plt.ioff()
    sns.set(style='white')
    vtypes = sorted(df["vtype"].unique(), reverse=True)
    if not callers:
        callers = sorted(df["caller"].unique())
    if not samples:
        samples = sorted(df["sample"].unique())
    if len(samples) >= len(callers):
        cats, groups = (samples, callers)
        data_dict = df.set_index(["sample", "caller", "vtype"]).T.to_dict()
    else:
        cats, groups = (callers, samples)
        data_dict = df.set_index(["caller", "sample", "vtype"]).T.to_dict()
    fig, axs = plt.subplots(len(vtypes) * len(groups), len(metrics))
    fig.text(.5, .95, title if title else "", horizontalalignment='center', size=14)
    for vi, vtype in enumerate(vtypes):
        sns.set_palette(sns.xkcd_palette([colors[vi]]))
        for gi, group in enumerate(groups):
            for mi, (metric, label) in enumerate(metrics):
                row_plots = axs if len(vtypes) * len(groups) == 1 else axs[vi * len(groups) + gi]
                cur_plot = row_plots if len(metrics) == 1 else row_plots[mi]
                vals, labels = [], []
                for cat in cats:
                    cur_data = data_dict.get((cat, group, vtype))
                    if cur_data:
                        vals.append(cur_data[metric])
                        labels.append(cur_data[label])
                cur_plot.barh(np.arange(len(vals)), vals)
                all_vals = []
                for k, d in data_dict.items():
                    if k[-1] == vtype:
                        for m in metrics:
                            all_vals.append(d[m[0]])
                metric_max = max(all_vals)
                cur_plot.set_xlim(0, metric_max)
                pad = 0.1 * metric_max
                for ai, (val, label) in enumerate(zip(vals, labels)):
                    cur_plot.annotate(label, (pad + (0 if max(vals) > metric_max / 2.0 else max(vals)),
                                              ai + 0.35), va='center', size=7)
                cur_plot.locator_params(nbins=len(cats) + 2, axis="y", tight=True)
                if mi == 0:
                    cur_plot.tick_params(axis='y', which='major', labelsize=8)
                    cur_plot.set_yticklabels(cats, size=8, va="bottom")
                    cur_plot.set_title("%s: %s" % (vtype, group), fontsize=12, loc="left")
                else:
                    cur_plot.get_yaxis().set_ticks([])
                if gi == len(groups) - 1:
                    cur_plot.tick_params(axis='x', which='major', labelsize=8)
                    cur_plot.get_xaxis().set_major_formatter(
                        mpl_ticker.FuncFormatter(lambda v, p: "%s%%" % (int(v) if round(v) == v else v)))
                    if vi == len(vtypes) - 1:
                        cur_plot.get_xaxis().set_label_text(metric_labels[metric], size=12)
                else:
                    cur_plot.get_xaxis().set_ticks([])
                    cur_plot.spines['bottom'].set_visible(False)
                cur_plot.spines['left'].set_visible(False)
                cur_plot.spines['top'].set_visible(False)
                cur_plot.spines['right'].set_visible(False)
    x, y = (6, len(vtypes) * len(groups) + 1 * 0.5 * len(cats)) if size is None else size
    fig.set_size_inches(x, y)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    plt.subplots_adjust(hspace=0.6)
    fig.savefig(out_file)

def create_from_csv(in_csv, config=None, outtype="png", title=None, size=None):
    df = pd.read_csv(in_csv)
    create(df, None, 0, config or {}, os.path.splitext(in_csv)[0], outtype, title,
           size)

def create(plot_data, header, ploti, sample_config, out_file_base, outtype="png",
           title=None, size=None):
    """Create plots of validation results for a sample, labeling prep strategies.
    """
    if mpl is None or plt is None or sns is None:
        not_found = ", ".join([x for x in ['mpl', 'plt', 'sns'] if eval(x) is None])
        logger.info("No validation plot. Missing imports: %s" % not_found)
        return None
    mpl.use('Agg', force=True)

    if header:
        df = pd.DataFrame(plot_data, columns=header)
    else:
        df = plot_data
    df["aligner"] = [get_aligner(x, sample_config) for x in df["sample"]]
    df["bamprep"] = [get_bamprep(x, sample_config) for x in df["sample"]]
    floors = get_group_floors(df, cat_labels)
    df["value.floor"] = [get_floor_value(x, cat, vartype, floors)
                         for (x, cat, vartype) in zip(df["value"], df["category"], df["variant.type"])]
    out = []
    for i, prep in enumerate(df["bamprep"].unique()):
        out.append(plot_prep_methods(df, prep, i + ploti, out_file_base, outtype, title, size))
    return out

cat_labels = {"concordant": "Concordant",
              "discordant-missing-total": "Discordant (missing)",
              "discordant-extra-total": "Discordant (extra)",
              "discordant-shared-total": "Discordant (shared)"}
vtype_labels = {"snp": "SNPs", "indel": "Indels"}
prep_labels = {}
caller_labels = {"ensemble": "Ensemble", "freebayes": "FreeBayes",
                 "gatk": "GATK Unified\nGenotyper", "gatk-haplotype": "GATK Haplotype\nCaller"}

def plot_prep_methods(df, prep, prepi, out_file_base, outtype, title=None,
                      size=None):
    """Plot comparison between BAM preparation methods.
    """
    samples = df[(df["bamprep"] == prep)]["sample"].unique()
    assert len(samples) >= 1, samples
    out_file = "%s-%s.%s" % (out_file_base, samples[0], outtype)
    df = df[df["category"].isin(cat_labels)]
    _seaborn(df, prep, prepi, out_file, title, size)
    return out_file

def _seaborn(df, prep, prepi, out_file, title=None, size=None):
    """Plot using seaborn wrapper around matplotlib.
    """
    plt.ioff()
    sns.set(style='dark')
    vtypes = df["variant.type"].unique()
    callers = sorted(df["caller"].unique())
    cats = _check_cats(["concordant", "discordant-missing-total",
                        "discordant-extra-total", "discordant-shared-total"],
                       vtypes, df, prep, callers)
    fig, axs = plt.subplots(len(vtypes), len(cats))
    width = 0.8
    for i, vtype in enumerate(vtypes):
        ax_row = axs[i] if len(vtypes) > 1 else axs
        for j, cat in enumerate(cats):
            vals, labels, maxval = _get_chart_info(df, vtype, cat, prep, callers)
            if len(cats) == 1:
                assert j == 0
                ax = ax_row
            else:
                ax = ax_row[j]
            if i == 0:
                ax.set_title(cat_labels[cat], size=14)
            ax.get_yaxis().set_ticks([])
            if j == 0:
                ax.set_ylabel(vtype_labels[vtype], size=14)
            ax.bar(np.arange(len(callers)), vals, width=width)
            ax.set_ylim(0, maxval)
            if i == len(vtypes) - 1:
                ax.set_xticks(np.arange(len(callers)) + width / 2.0)
                ax.set_xticklabels([caller_labels.get(x, x).replace("__", "\n") if x else ""
                                    for x in callers], size=8, rotation=45)
            else:
                ax.get_xaxis().set_ticks([])
            _annotate(ax, labels, vals, np.arange(len(callers)), width)
    fig.text(.5, .95, prep_labels.get(prep, "") if title is None else title, horizontalalignment='center', size=16)
    fig.subplots_adjust(left=0.05, right=0.95, top=0.87, bottom=0.15, wspace=0.1, hspace=0.1)
    x, y = (10, 5) if size is None else size
    fig.set_size_inches(x, y)
    fig.savefig(out_file)

def _check_cats(cats, vtypes, df, prep, callers):
    """Only include categories in the final output if they have values.
    """
    out = []
    for cat in cats:
        all_vals = []
        for vtype in vtypes:
            vals, labels, maxval = _get_chart_info(df, vtype, cat, prep, callers)
            all_vals.extend(vals)
        if sum(all_vals) / float(len(all_vals)) > 2:
            out.append(cat)
    if len(out) == 0:
        return cats
    else:
        return out

def _get_chart_info(df, vtype, cat, prep, callers):
    """Retrieve values for a specific variant type, category and prep method.
    """
    maxval_raw = max(list(df["value.floor"]))
    curdf = df[(df["variant.type"] == vtype) & (df["category"] == cat)
               & (df["bamprep"] == prep)]
    vals = []
    labels = []
    for c in callers:
        row = curdf[df["caller"] == c]
        if len(row) > 0:
            vals.append(list(row["value.floor"])[0])
            labels.append(list(row["value"])[0])
        else:
            vals.append(1)
            labels.append("")
    return vals, labels, maxval_raw

def _annotate(ax, annotate, height, left, width):
    """Annotate axis with labels.
    """
    annotate_yrange_factor = 0.010
    xticks = np.array(left) + width / 2.0
    ymin, ymax = ax.get_ylim()
    yrange = ymax - ymin

    # Reset ymax and ymin so there's enough room to see the annotation of
    # the top-most
    if ymax > 0:
        ymax += yrange * 0.15
    if ymin < 0:
        ymin -= yrange * 0.15
    ax.set_ylim(ymin, ymax)
    yrange = ymax - ymin

    offset_ = yrange * annotate_yrange_factor
    if isinstance(annotate, collections.Iterable):
        annotations = map(str, annotate)
    else:
        annotations = ['%.3f' % h if type(h) is np.float_ else str(h)
                       for h in height]
    for x, h, annotation in zip(xticks, height, annotations):
        # Adjust the offset to account for negative bars
        offset = offset_ if h >= 0 else -1 * offset_
        verticalalignment = 'bottom' if h >= 0 else 'top'

        if len(str(annotation)) > 6:
            size = 7
        elif len(str(annotation)) > 5:
            size = 8
        else:
            size = 10
        # Finally, add the text to the axes
        ax.annotate(annotation, (x, h + offset),
                    verticalalignment=verticalalignment,
                    horizontalalignment='center',
                    size=size)

def _ggplot(df, out_file):
    """Plot faceted items with ggplot wrapper on top of matplotlib.
    XXX Not yet functional
    """
    import ggplot as gg
    df["variant.type"] = [vtype_labels[x] for x in df["variant.type"]]
    df["category"] = [cat_labels[x] for x in df["category"]]
    df["caller"] = [caller_labels.get(x, None) for x in df["caller"]]
    p = (gg.ggplot(df, gg.aes(x="caller", y="value.floor")) + gg.geom_bar()
         + gg.facet_wrap("variant.type", "category")
         + gg.theme_seaborn())
    gg.ggsave(p, out_file)

def get_floor_value(x, cat, vartype, floors):
    """Modify values so all have the same relative scale for differences.

    Using the chosen base heights, adjusts an individual sub-plot to be consistent
    relative to that height.
    """
    all_base = floors[vartype]
    cur_max = floors[(cat, vartype)]
    if cur_max > all_base:
        diff = cur_max - all_base
        x = max(1, x - diff)
    return x

def get_group_floors(df, cat_labels):
    """Retrieve the floor for a given row of comparisons, creating a normalized set of differences.

    We need to set non-zero floors so large numbers (like concordance) don't drown out small
    numbers (like discordance). This defines the height for a row of comparisons as either
    the minimum height of any sub-plot, or the maximum difference between higher and lower
    (plus 10%).
    """
    group_maxes = collections.defaultdict(list)
    group_diffs = collections.defaultdict(list)
    diff_pad = 0.1  # 10% padding onto difference to avoid large numbers looking like zero
    for name, group in df.groupby(["category", "variant.type"]):
        label, stype = name
        if label in cat_labels:
            diff = max(group["value"]) - min(group["value"])
            group_diffs[stype].append(diff + int(diff_pad * diff))
            group_maxes[stype].append(max(group["value"]))
        group_maxes[name].append(max(group["value"]))
    out = {}
    for k, vs in group_maxes.iteritems():
        if k in group_diffs:
            out[k] = max(max(group_diffs[stype]), min(vs))
        else:
            out[k] = min(vs)
    return out

def get_aligner(x, config):
    return utils.get_in(config, ("algorithm", "aligner"), "")

def get_bamprep(x, config):
    params = bamprep._get_prep_params({"config": {"algorithm": config.get("algorithm", {})}})
    if params["realign"] == "gatk" and params["recal"] == "gatk":
        return "gatk"
    elif not params["realign"] and not params["recal"]:
        return "none"
    elif not params.get("recal") or not params.get("realign"):
        return "mixed"
    else:
        return ""

# ## Frequency plots

def facet_freq_plot(freq_csv, caller):
    """Prepare a facet plot of frequencies stratified by variant type and status (TP, FP, FN).

    Makes a nice plot with the output from validate.freq_summary
    """
    out_file = "%s.png" % os.path.splitext(freq_csv)[0]
    plt.ioff()
    sns.set(style='dark')
    df = pd.read_csv(freq_csv)
    g = sns.FacetGrid(df, row="vtype", col="valclass", margin_titles=True,
                      col_order=["TP", "FN", "FP"], row_order=["snp", "indel"],
                      sharey=False)
    g.map(plt.hist, "freq", bins=20, align="left")
    g.set(xlim=(0.0, 1.0))
    g.fig.set_size_inches(8, 6)
    g.fig.text(.05, .97, caller, horizontalalignment='center', size=14)
    g.fig.savefig(out_file)
