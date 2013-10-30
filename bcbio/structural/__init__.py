"""Detect structural variation in genomes using high-throughput sequencing data.
"""
import collections
import copy

from bcbio.structural import lumpy

_CALLERS = {"lumpy": lumpy.run}

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

def _combine_multiple_svcallers(samples):
    """
    """
    by_bam = collections.defaultdict(list)
    for x in samples:
        by_bam[x[0]["work_bam"]].append(x[0])
    out = []
    for grouped_calls in by_bam.itervalues():
        def orig_svcaller_order(x):
            return _get_svcallers(x).index(x["config"]["algorithm"]["svcaller_active"])
        sorted_svcalls = sorted([x for x in grouped_calls if "sv" in x],
                                key=orig_svcaller_order)
        final_calls = [x["sv"] for x in sorted_svcalls]
        final = grouped_calls[0]
        final["sv"] = final_calls
        del final["config"]["algorithm"]["svcaller_active"]
        out.append([final])
    return out

def run(samples, run_parallel):
    """Run structural variation detection using configured methods.
    """
    to_process = []
    extras = []
    for data in (xs[0] for xs in samples):
        ready_data = _handle_multiple_svcallers(data)
        if len(ready_data) > 0:
            for x in ready_data:
                to_process.append([x])
        else:
            extras.append([data])
    processed = run_parallel("detect_sv", to_process)
    return extras + _combine_multiple_svcallers(processed)

def detect_sv(data):
    """Top level parallel target for examining structural variation.
    """
    svcaller = data["config"]["algorithm"].get("svcaller_active")
    if svcaller:
        data["sv"] = _CALLERS[svcaller](data)
    return [[data]]
