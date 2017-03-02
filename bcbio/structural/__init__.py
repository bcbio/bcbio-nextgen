"""Detect structural variation in genomes using high-throughput sequencing data.
"""
import collections
import copy
import operator

import toolz as tz

from bcbio import utils
from bcbio.cwl import cwlutils
from bcbio.pipeline import datadict as dd
from bcbio.structural import (battenberg, cn_mops, cnvkit, delly, gridss,
                              lumpy, manta, metasv, prioritize, plot,
                              seq2c, validate, wham)
from bcbio.variation import vcfutils

# Stratify callers by stage -- see `run` documentation below for definitions
_CALLERS = {
  "precall": {"seq2c": seq2c.precall},
  "initial": {"cnvkit": cnvkit.run},
  "standard": {"cn.mops": cn_mops.run, "manta": manta.run, "cnvkit": cnvkit.run,
               "delly": delly.run, "lumpy": lumpy.run, "wham": wham.run,
               "battenberg": battenberg.run, "seq2c": seq2c.run, "gridss": gridss.run},
  "ensemble": {"metasv": metasv.run,
               "prioritize": prioritize.run}}
_NEEDS_BACKGROUND = set(["cn.mops"])
_GLOBAL_BATCHING = set(["seq2c"])

def _get_callers(items, stage):
    """Retrieve available callers for the provided stage.

    Handles special cases like CNVkit that can be in initial or standard
    depending on if fed into Lumpy analysis.
    """
    callers = _CALLERS[stage]
    if "cnvkit" in callers:
        has_lumpy = any("lumpy" in get_svcallers(d) or lumpy in d["config"]["algorithm"].get("svcaller_orig", [])
                        for d in items)
        if has_lumpy and any("lumpy_usecnv" in dd.get_tools_on(d) for d in items):
            if stage != "initial":
                del callers["cnvkit"]
        else:
            if stage != "standard":
                del callers["cnvkit"]
    return callers

def get_svcallers(data):
    svs = data["config"]["algorithm"].get("svcaller")
    if svs is None:
        svs = []
    elif isinstance(svs, basestring):
        svs = [svs]
    return svs

def _handle_multiple_svcallers(data, stage):
    """Retrieve configured structural variation caller, handling multiple.
    """
    svs = get_svcallers(data)
    # special cases -- prioritization
    if stage == "ensemble" and dd.get_svprioritize(data):
        svs.append("prioritize")
    out = []
    for svcaller in svs:
        if svcaller in _get_callers([data], stage):
            base = copy.deepcopy(data)
            base["config"]["algorithm"]["svcaller"] = svcaller
            base["config"]["algorithm"]["svcaller_orig"] = svs
            out.append(base)
    return out

def finalize_sv(samples, config):
    """Combine results from multiple sv callers into a single ordered 'sv' key.
    """
    by_bam = collections.OrderedDict()
    for x in samples:
        batch = dd.get_batch(x) or [dd.get_sample_name(x)]
        try:
            by_bam[x["align_bam"], tuple(batch)].append(x)
        except KeyError:
            by_bam[x["align_bam"], tuple(batch)] = [x]
    by_batch = collections.OrderedDict()
    lead_batches = {}
    for grouped_calls in by_bam.values():
        def orig_svcaller_order(x):
            orig_callers = tz.get_in(["config", "algorithm", "svcaller_orig"], x)
            cur_caller = tz.get_in(["config", "algorithm", "svcaller"], x)
            return orig_callers.index(cur_caller)
        sorted_svcalls = sorted([x for x in grouped_calls if "sv" in x],
                                key=orig_svcaller_order)
        final = grouped_calls[0]
        if len(sorted_svcalls) > 0:
            final["sv"] = reduce(operator.add, [x["sv"] for x in sorted_svcalls])
        final["config"]["algorithm"]["svcaller"] = final["config"]["algorithm"].pop("svcaller_orig")
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

def batch_for_sv(samples):
    """Prepare a set of samples for parallel structural variant calling.

    CWL input target -- groups samples into batches and structural variant
    callers for parallel processing.
    """
    to_process, extras, background = _batch_split_by_sv(samples, "standard")
    out = [cwlutils.samples_to_records(xs) for xs in to_process.values()] + extras
    return out

def _batch_split_by_sv(samples, stage):
    to_process = collections.OrderedDict()
    extras = []
    background = []
    for data in (utils.to_single_data(x) for x in samples):
        ready_data = _handle_multiple_svcallers(data, stage)
        if len(ready_data) > 0:
            background.append(data)
            for x in ready_data:
                svcaller = tz.get_in(["config", "algorithm", "svcaller"], x)
                batch = dd.get_batch(x) or dd.get_sample_name(x)
                if stage in ["precall", "ensemble"]:  # no batching for precall or ensemble methods
                    batch = "%s-%s" % (dd.get_sample_name(x), batch)
                elif svcaller in _GLOBAL_BATCHING:  # All samples batched together for analyses
                    batch = "all"
                batches = batch if isinstance(batch, (list, tuple)) else [batch]
                for b in batches:
                    try:
                        to_process[(svcaller, b)].append(x)
                    except KeyError:
                        to_process[(svcaller, b)] = [x]
        else:
            extras.append([data])
    return to_process, extras, background

def run(samples, run_parallel, stage):
    """Run structural variation detection.

    The stage indicates which level of structural variant calling to run.
      - precall, perform initial sample based assessment of samples
      - initial, callers that can be used in subsequent structural variation steps (cnvkit -> lumpy)
      - standard, regular batch calling
      - ensemble, post-calling, combine other callers or prioritize results
    """
    to_process, extras, background = _batch_split_by_sv(samples, stage)
    processed = run_parallel("detect_sv", ([xs, background, stage]
                                           for xs in to_process.values()))
    finalized = (run_parallel("finalize_sv", [([xs[0] for xs in processed], processed[0][0]["config"])])
                 if len(processed) > 0 else [])
    return extras + finalized

def detect_sv(items, all_items=None, stage="standard"):
    """Top level parallel target for examining structural variation.
    """
    svcaller = items[0]["config"]["algorithm"].get("svcaller")
    caller_fn = _get_callers(items, stage).get(svcaller)
    out = []
    if svcaller and caller_fn:
        if (all_items and svcaller in _NEEDS_BACKGROUND and
                not vcfutils.is_paired_analysis([x.get("align_bam") for x in items], items)):
            names = set([dd.get_sample_name(x) for x in items])
            background = [x for x in all_items if dd.get_sample_name(x) not in names]
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
    """Use more resources (up to available limits) if we have multiple QC samples/svcallers.
    """
    machines = []
    for data in (xs[0] for xs in items):
        machines.append(max(1, len(get_svcallers(data)), len(dd.get_algorithm_qc(data))))
    return sum(machines)
