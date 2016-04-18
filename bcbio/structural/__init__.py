"""Detect structural variation in genomes using high-throughput sequencing data.
"""
import collections
import copy
import operator

import toolz as tz

from bcbio.pipeline import datadict as dd
from bcbio.structural import (battenberg, cn_mops, cnvkit, delly,
                              lumpy, manta, metasv, prioritize, plot, validate, wham)
from bcbio.variation import vcfutils

# Stratify callers by stage -- see `run` documentation below for definitions
_CALLERS = {
  "standard": {"cn.mops": cn_mops.run, "manta": manta.run,
               "delly": delly.run, "lumpy": lumpy.run, "wham": wham.run,
               "cnvkit": cnvkit.run, "battenberg": battenberg.run},
  "ensemble": {"metasv": metasv.run,
               "prioritize": prioritize.run}}
_NEEDS_BACKGROUND = set(["cn.mops"])

def _get_svcallers(data):
    svs = data["config"]["algorithm"].get("svcaller")
    if svs is None:
        svs = []
    elif isinstance(svs, basestring):
        svs = [svs]
    return svs

def _handle_multiple_svcallers(data, stage):
    """Retrieve configured structural variation caller, handling multiple.
    """
    svs = _get_svcallers(data)
    # special cases -- prioritization
    if stage == "ensemble" and dd.get_svprioritize(data):
        svs.append("prioritize")
    out = []
    for svcaller in svs:
        if svcaller in _CALLERS[stage]:
            base = copy.deepcopy(data)
            base["config"]["algorithm"]["svcaller"] = svs
            base["config"]["algorithm"]["svcaller_active"] = svcaller
            out.append(base)
    return out

def finalize_sv(samples, config):
    """Combine results from multiple sv callers into a single ordered 'sv' key.
    """
    by_bam = collections.OrderedDict()
    for x in samples:
        try:
            by_bam[x["align_bam"]].append(x)
        except KeyError:
            by_bam[x["align_bam"]] = [x]
    by_batch = collections.OrderedDict()
    lead_batches = {}
    for grouped_calls in by_bam.values():
        def orig_svcaller_order(x):
            return _get_svcallers(x).index(x["config"]["algorithm"]["svcaller_active"])
        sorted_svcalls = sorted([x for x in grouped_calls if "sv" in x],
                                key=orig_svcaller_order)
        final = grouped_calls[0]
        if len(sorted_svcalls) > 0:
            final["sv"] = reduce(operator.add, [x["sv"] for x in sorted_svcalls])
        del final["config"]["algorithm"]["svcaller_active"]
        batch = dd.get_batch(final) or dd.get_sample_name(final)
        batches = batch if isinstance(batch, (list, tuple)) else [batch]
        lead_batches[dd.get_sample_name(final)] = batches[0]
        for batch in batches:
            try:
                by_batch[batch].append(final)
            except KeyError:
                by_batch[batch] = [final]
    out = []
    for batch, items in by_batch.items():
        if any("svplots" in dd.get_tools_on(d) for d in items):
            plot_items = plot.by_regions(items)
        else:
            plot_items = items
        for data in plot_items:
            if lead_batches[dd.get_sample_name(data)] == batch:
                out.append([data])
    return out

def validate_sv(data):
    """Validate structural variant calls for a sample.
    """
    return [[validate.evaluate(data)]]

def run(samples, run_parallel, stage):
    """Run structural variation detection.

    The stage indicates which level of structural variant calling to run.
      - standard, regular batch calling
      - ensemble, post-calling, combine other callers or prioritize results
    """
    to_process = collections.OrderedDict()
    extras = []
    background = []
    for data in (xs[0] for xs in samples):
        ready_data = _handle_multiple_svcallers(data, stage)
        if len(ready_data) > 0:
            background.append(data)
            for x in ready_data:
                svcaller = x["config"]["algorithm"].get("svcaller_active")
                batch = dd.get_batch(x) or dd.get_sample_name(x)
                if stage == "ensemble":  # no batching for ensemble methods
                    batch = "%s-%s" % (dd.get_sample_name(x), batch)
                batches = batch if isinstance(batch, (list, tuple)) else [batch]
                for b in batches:
                    try:
                        to_process[(svcaller, b)].append(x)
                    except KeyError:
                        to_process[(svcaller, b)] = [x]
        else:
            extras.append([data])
    processed = run_parallel("detect_sv", ([xs, background, xs[0]["config"], stage]
                                           for xs in to_process.values()))
    finalized = (run_parallel("finalize_sv", [([xs[0] for xs in processed], processed[0][0]["config"])])
                 if len(processed) > 0 else [])
    return extras + finalized

def detect_sv(items, all_items, config, stage):
    """Top level parallel target for examining structural variation.
    """
    svcaller = config["algorithm"].get("svcaller_active")
    caller_fn = _CALLERS[stage].get(svcaller)
    out = []
    if svcaller and caller_fn:
        if (svcaller in _NEEDS_BACKGROUND and
                not vcfutils.is_paired_analysis([x.get("align_bam") for x in items], items)):
            names = set([tz.get_in(["rgnames", "sample"], x) for x in items])
            background = [x for x in all_items if tz.get_in(["rgnames", "sample"], x) not in names]
            for svdata in caller_fn(items, background):
                out.append([svdata])
        else:
            for svdata in caller_fn(items):
                out.append([svdata])
    else:
        for data in items:
            out.append([data])
    return out

# ## configuration

def parallel_multiplier(items):
    """Use more resources (up to available limits) if we have multiple SV callers.
    """
    return max([1] + [len(_get_svcallers(xs[0])) for xs in items])
