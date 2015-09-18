"""Detect structural variation in genomes using high-throughput sequencing data.
"""
import collections
import copy
import operator

import toolz as tz

from bcbio.pipeline import datadict as dd
from bcbio.structural import (battenberg, cn_mops, cnvkit, delly, ensemble,
                              lumpy, manta, metasv, plot, validate, wham)
from bcbio.variation import vcfutils

_CALLERS = {}
_SOMATIC_CALLERS = {"manta": manta.run}
_BATCH_CALLERS = {"cn.mops": cn_mops.run, "cnvkit": cnvkit.run,
                  "delly": delly.run, "lumpy": lumpy.run, "wham": wham.run,
                  "battenberg": battenberg.run}
_ENSEMBLE_CALLERS = {"metasv": metasv.run}
_NEEDS_BACKGROUND = set(["cn.mops"])
_INITIAL_CALLERS = set(["battenberg", "cnvkit"])

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
        if svcaller not in _ENSEMBLE_CALLERS:
            base = copy.deepcopy(data)
            base["config"]["algorithm"]["svcaller_active"] = svcaller
            out.append(base)
    return out

def finalize_sv(samples, config, initial_only=False):
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
            if not initial_only:
                for caller in (c for c in _get_svcallers(final) if c in _ENSEMBLE_CALLERS):
                    final_calls = _ENSEMBLE_CALLERS[caller](final_calls, final)
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

def run(samples, run_parallel, initial_only=False):
    """Run structural variation detection.

    initial_only indicates we should run structural variation inputs, like
    CNV calling, we can use to inform low frequency variant calling.
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
                # reset SV information if we're running a second pass SV call
                if "sv" in x:
                    del x["sv"]
                batch = dd.get_batch(x)
                paired = vcfutils.get_paired_phenotype(x)
                if ((svcaller in _BATCH_CALLERS and batch) or
                      (svcaller in _SOMATIC_CALLERS and paired and batch)):
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
    processed = run_parallel("detect_sv", ([xs, background, xs[0]["config"], initial_only]
                                           for xs in to_process.values()))
    finalized = (run_parallel("finalize_sv", [([xs[0] for xs in processed], processed[0][0]["config"],
                                               initial_only)])
                 if len(processed) > 0 else [])
    return extras + finalized

def detect_sv(items, all_items, config, initial_only=False):
    """Top level parallel target for examining structural variation.
    """
    svcaller = config["algorithm"].get("svcaller_active")
    out = []
    if svcaller and (not initial_only or svcaller in _INITIAL_CALLERS):
        if svcaller in _CALLERS:
            assert len(items) == 1
            data = items[0]
            data["sv"] = _CALLERS[svcaller](data)
            out.append([data])
        elif svcaller in _BATCH_CALLERS or svcaller in _SOMATIC_CALLERS:
            if (svcaller in _NEEDS_BACKGROUND and
                  not vcfutils.is_paired_analysis([x.get("align_bam") for x in items], items)):
                names = set([tz.get_in(["rgnames", "sample"], x) for x in items])
                background = [x for x in all_items if tz.get_in(["rgnames", "sample"], x) not in names]
                for svdata in _BATCH_CALLERS[svcaller](items, background):
                    out.append([svdata])
            else:
                for svdata in _BATCH_CALLERS.get(svcaller, _SOMATIC_CALLERS.get(svcaller))(items):
                    out.append([svdata])
        else:
            raise ValueError("Unexpected structural variant caller: %s" % svcaller)
    else:
        for data in items:
            out.append([data])
    return out

# ## configuration

def parallel_multiplier(items):
    """Use more resources (up to available limits) if we have multiple SV callers.
    """
    return max([1] + [len(_get_svcallers(xs[0])) for xs in items])
