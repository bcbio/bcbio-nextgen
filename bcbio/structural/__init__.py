"""Detect structural variation in genomes using high-throughput sequencing data.
"""
import collections
import copy
import operator

import toolz as tz

from bcbio.pipeline import datadict as dd
from bcbio.structural import (battenberg, cn_mops, cnvkit, delly, ensemble,
                              lumpy, plot, validate, wham)
from bcbio.variation import vcfutils

_CALLERS = {}
_BATCH_CALLERS = {"cn.mops": cn_mops.run, "cnvkit": cnvkit.run,
                  "delly": delly.run, "lumpy": lumpy.run, "wham": wham.run,
                  "battenberg": battenberg.run}
_NEEDS_BACKGROUND = set(["cn.mops"])

def _get_svcallers(data):
    svs = data["config"]["algorithm"].get("svcaller")
    if svs is None:
        svs = []
    elif isinstance(svs, basestring):
        svs = [svs]
    return svs

def _handle_multiple_svcallers(data):
    """Retrieve configured structural variation caller, handling multiple.
    """
    svs = _get_svcallers(data)
    out = []
    for svcaller in svs:
        base = copy.deepcopy(data)
        base["config"]["algorithm"]["svcaller_active"] = svcaller
        out.append(base)
    return out

def finalize_sv(samples, config):
    """Combine results from multiple sv callers into a single ordered 'sv' key.

    Handles ensemble calling and plotting of results.
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
            final_calls = reduce(operator.add, [x["sv"] for x in sorted_svcalls])
            final_calls = ensemble.summarize(final_calls, final, grouped_calls)
            final_calls = validate.evaluate(final, final_calls)
            final["sv"] = final_calls
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

def run(samples, run_parallel):
    """Run structural variation detection using configured methods.
    """
    to_process = collections.OrderedDict()
    extras = []
    background = []
    for data in (xs[0] for xs in samples):
        ready_data = _handle_multiple_svcallers(data)
        if len(ready_data) > 0:
            background.append(data)
            for x in ready_data:
                svcaller = x["config"]["algorithm"].get("svcaller_active")
                batch = dd.get_batch(x)
                if svcaller in _BATCH_CALLERS and batch:
                    batches = batch if isinstance(batch, (list, tuple)) else [batch]
                    for b in batches:
                        try:
                            to_process[(svcaller, b)].append(x)
                        except KeyError:
                            to_process[(svcaller, b)] = [x]
                else:
                    to_process[(svcaller, dd.get_sample_name(x))] = [x]
        else:
            extras.append([data])
    processed = run_parallel("detect_sv", ([xs, background, xs[0]["config"]] for xs in to_process.values()))
    finalized = (run_parallel("finalize_sv", [([xs[0] for xs in processed], processed[0][0]["config"])])
                 if len(processed) > 0 else [])
    return extras + finalized

def detect_sv(items, all_items, config):
    """Top level parallel target for examining structural variation.
    """
    svcaller = config["algorithm"].get("svcaller_active")
    out = []
    if svcaller:
        if svcaller in _CALLERS:
            assert len(items) == 1
            data = items[0]
            data["sv"] = _CALLERS[svcaller](data)
            out.append([data])
        elif svcaller in _BATCH_CALLERS:
            if (svcaller in _NEEDS_BACKGROUND and
                  not vcfutils.is_paired_analysis([x.get("align_bam") for x in items], items)):
                names = set([tz.get_in(["rgnames", "sample"], x) for x in items])
                background = [x for x in all_items if tz.get_in(["rgnames", "sample"], x) not in names]
                for svdata in _BATCH_CALLERS[svcaller](items, background):
                    out.append([svdata])
            else:
                for svdata in _BATCH_CALLERS[svcaller](items):
                    out.append([svdata])
        else:
            raise ValueError("Unexpected structural variant caller: %s" % svcaller)
    else:
        out.append(items)
    return out
