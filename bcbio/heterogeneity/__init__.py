"""Calculate heterogeneity and sub-clonal populations for complex input samples.

Use allele frequencies, copy number calls and structural variants to infer
sub-clonal populations within a potentially mixed population. This attempts
to infer these sub-clones to help improve variant calls and interpretation
especially in complex cancer samples.
"""
import collections

from bcbio.heterogeneity import bubbletree, phylowgs, theta
from bcbio.pipeline import datadict as dd
from bcbio.variation import vcfutils

def _get_calls(data, cnv_only=False):
    """Retrieve calls, organized by name, to use for heterogeneity analysis.
    """
    cnvs_supported = set(["cnvkit", "battenberg"])
    out = {}
    for sv in data.get("sv", []):
        if not cnv_only or sv["variantcaller"] in cnvs_supported:
            out[sv["variantcaller"]] = sv
    return out

def _get_variants(data):
    """Retrieve set of variant calls to use for heterogeneity analysis.
    """
    supported = ["vardict", "vardict-java", "vardict-perl", "freebayes", "mutect"]
    out = []
    for v in data.get("variants", []):
        if v["variantcaller"] in supported:
            out.append((supported.index(v["variantcaller"]), v))
    out.sort()
    return [xs[1] for xs in out]

def _ready_for_het_analysis(items):
    """Check if a sample has input information for heterogeneity analysis.

    We currently require a tumor/normal sample containing both CNV and variant calls.
    """
    paired = vcfutils.get_paired_bams([dd.get_align_bam(d) for d in items], items)
    has_het = any(dd.get_hetcaller(d) for d in items)
    if has_het and paired:
        return _get_variants(paired.tumor_data) and _get_calls(paired.tumor_data, cnv_only=True)

def _get_batches(data):
    batches = dd.get_batch(data) or dd.get_sample_name(data)
    if not isinstance(batches, (list, tuple)):
        batches = [batches]
    return batches

def _group_by_batches(items):
    out = collections.OrderedDict()
    for data in (xs[0] for xs in items):
        for b in _get_batches(data):
            try:
                out[b].append(data)
            except KeyError:
                out[b] = [data]
    return out

def _get_hetcallers(items):
    out = set([])
    for d in items:
        hetcaller = dd.get_hetcaller(d)
        if hetcaller:
            out = out.union(set(hetcaller))
    return sorted(list(out))

def estimate(items, batch, config):
    """Estimate heterogeneity for a pair of tumor/normal samples. Run in parallel.
    """
    hetcallers = {"theta": theta.run,
                  "phylowgs": phylowgs.run,
                  "bubbletree": bubbletree.run}
    paired = vcfutils.get_paired_bams([dd.get_align_bam(d) for d in items], items)
    calls = _get_calls(paired.tumor_data)
    variants = _get_variants(paired.tumor_data)
    het_info = []
    for hetcaller in _get_hetcallers(items):
        try:
            hetfn = hetcallers[hetcaller]
        except KeyError:
            hetfn = None
            print "%s not yet implemented" % hetcaller
        if hetfn:
            hetout = hetfn(variants[0], calls, paired)
            if hetout:
                het_info.append(hetout)
    out = []
    for data in items:
        if batch == _get_batches(data)[0]:
            if dd.get_sample_name(data) == paired.tumor_name:
                if het_info:
                    data["heterogeneity"] = het_info
            out.append([data])
    return out

def run(items, run_parallel):
    """Top level entry point for calculating heterogeneity, handles organization and job distribution.
    """
    to_process = []
    extras = []
    for batch, cur_items in _group_by_batches(items).items():
        if _ready_for_het_analysis(cur_items):
            to_process.append((batch, cur_items))
        else:
            for data in cur_items:
                extras.append([data])
    processed = run_parallel("heterogeneity_estimate", ([xs, b, xs[0]["config"]] for b, xs in to_process))
    return extras + processed
