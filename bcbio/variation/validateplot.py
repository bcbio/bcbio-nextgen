"""Plot validation results from variant calling comparisons.

Handles data normalization and plotting, emphasizing comparisons on methodology
differences.
"""
import bisect
import collections
import math

import numpy as np
try:
    import prettyplotlib as ppl
    import pandas as pd
except ImportError:
    gg, pd, ppl = None, None, None

from bcbio.variation import bamprep

def create(plot_data, header, ploti, sample_config, out_file_base):
    """Create plots of validation results for a sample, labeling prep strategies.
    """
    if pd is None or ppl is None:
        return None
    df = pd.DataFrame(plot_data, columns=header)
    df["aligner"] = [get_aligner(x, sample_config) for x in df["sample"]]
    df["bamprep"] = [get_bamprep(x, sample_config) for x in df["sample"]]
    floors = get_group_floors(df)
    df["value.floor"] = [get_floor_value(x, cat, vartype, floors)
                         for (x, cat, vartype) in zip(df["value"], df["category"], df["variant.type"])]
    out = []
    for i, prep in enumerate(df["bamprep"].unique()):
        out.append(plot_prep_methods(df, prep, i + ploti, out_file_base))
    return out

cat_labels = {"concordant": "Concordant",
              "discordant-missing-total": "Discordant (missing)",
              "discordant-extra-total": "Discordant (extra)",
              "discordant-shared-total": "Discordant (shared)"}
vtype_labels = {"snp": "SNPs", "indel": "Indels"}
prep_labels = {"gatk": "GATK best-practice BAM preparation (recalibration, realignment)",
               "none": "Minimal BAM preparation (samtools de-duplication only)"}
caller_labels = {"ensemble": "Ensemble", "freebayes": "FreeBayes",
                 "gatk": "GATK Unified\nGenotyper", "gatk-haplotype": "GATK Haplotype\nCaller"}

def plot_prep_methods(df, prep, prepi, out_file_base):
    """Plot comparison between BAM preparation methods.
    """
    samples = df[(df["bamprep"] == prep)]["sample"].unique()
    assert len(samples) >= 1, samples
    out_file = "%s-%s.pdf" % (out_file_base, samples[0])
    df = df[df["category"].isin(cat_labels)]
    _prettyplot(df, prep, prepi, out_file)
    return out_file

def _prettyplot(df, prep, prepi, out_file):
    """Plot using prettyplot wrapper around matplotlib.
    """
    cats = ["concordant", "discordant-missing-total",
            "discordant-extra-total", "discordant-shared-total"]
    vtypes = df["variant.type"].unique()
    fig, axs = ppl.subplots(len(vtypes), len(cats))
    callers = sorted(df["caller"].unique())
    width = 0.8
    for i, vtype in enumerate(vtypes):
        for j, cat in enumerate(cats):
            ax = axs[i][j]
            if i == 0:
                ax.set_title(cat_labels[cat], size=14)
            ax.get_yaxis().set_ticks([])
            if j == 0:
                ax.set_ylabel(vtype_labels[vtype], size=14)
            vals, labels, maxval = _get_chart_info(df, vtype, cat, prep, callers)
            ppl.bar(ax, np.arange(len(callers)), vals,
                    color=ppl.colors.set2[prepi], width=width)
            ax.set_ylim(0, maxval)
            if i == len(vtypes) - 1:
                ax.set_xticks(np.arange(len(callers)) + width / 2.0)
                ax.set_xticklabels([caller_labels.get(x, x) for x in callers], size=8, rotation=45)
            else:
                ax.get_xaxis().set_ticks([])
            _annotate(ax, labels, vals, np.arange(len(callers)), width)
    fig.text(.5, .95, prep_labels[prep], horizontalalignment='center', size=16)
    fig.subplots_adjust(left=0.05, right=0.95, top=0.87, bottom=0.15, wspace=0.1, hspace=0.1)
    #fig.tight_layout()
    fig.set_size_inches(10, 5)
    fig.savefig(out_file)

def _get_chart_info(df, vtype, cat, prep, callers):
    """Retrieve values for a specific variant type, category and prep method.
    """
    maxval_raw = max(list(df["value.floor"]))
    norm_ylim = 1000.0  # ceil to make plots more comparable
    maxval = math.ceil(maxval_raw / norm_ylim) * norm_ylim
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
    return vals, labels, maxval

def _annotate(ax, annotate, height, left, width):
    """Annotate axis with labels. Adjusted from prettyplotlib to be more configurable.
    Needed to adjust label size.
    """
    annotate_yrange_factor = 0.025
    xticks = np.array(left) + width / 2.0
    ymin, ymax = ax.get_ylim()
    yrange = ymax - ymin

    # Reset ymax and ymin so there's enough room to see the annotation of
    # the top-most
    if ymax > 0:
        ymax += yrange * 0.1
    if ymin < 0:
        ymin -= yrange * 0.1
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

        # Finally, add the text to the axes
        ax.annotate(annotation, (x, h + offset),
                    verticalalignment=verticalalignment,
                    horizontalalignment='center',
                    size=10,
                    color=ppl.colors.almost_black)

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
    base = floors[(cat, vartype)]
    #print cat, vartype, x, base
    return x - base

def get_group_floors(df):
    """Floor values to nearest 5,000 for each category.
    """
    fudge = 500
    floors = {}
    floor_vals = [x * 5e3 for x in range(50)]
    for name, group in df.groupby(["category", "variant.type"]):
        min_value = max(0, min(group["value"]) - fudge)
        floors[name] = int(floor_vals[bisect.bisect(floor_vals, min_value) - 1])
    return floors

def get_aligner(x, config):
    return config["algorithm"].get("aligner", "")

def get_bamprep(x, config):
    params = bamprep._get_prep_params({"config": {"algorithm": config["algorithm"]}})
    if params["realign"] == "gatk" and params["recal"] == "gatk":
        return "gatk"
    elif not params["realign"] and not params["recal"]:
        return "none"
    else:
        raise ValueError("Unexpected bamprep approach: %s" % params)
