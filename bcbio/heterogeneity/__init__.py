"""Calculate heterogeneity and sub-clonal populations for complex input samples.

Use allele frequencies, copy number calls and structural variants to infer
sub-clonal populations within a potentially mixed population. This attempts
to infer these sub-clones to help improve variant calls and interpretation
especially in complex cancer samples.
"""
import collections

from bcbio.heterogeneity import theta
from bcbio.pipeline import datadict as dd
from bcbio.variation import vcfutils

def _get_cnvs(data):
    """Retrieve CNV calls to use for heterogeneity analysis.
    """
    supported = set(["cnvkit"])
    out = []
    for sv in data.get("sv", []):
        if sv["variantcaller"] in supported:
            out.append(sv)
    return out

def _get_variants(data):
    """Retrieve set of variant calls to use for heterogeneity analysis.
    """
    supported = set(["ensemble", "vardict", "freebayes", "mutect"])
    out = []
    for v in data.get("variants", []):
        if v["variantcaller"] in supported:
            out.append(v)
    return out

def _ready_for_het_analysis(items):
    """Check if a sample has input information for heterogeneity analysis.

    We currently require a tumor/normal sample containing both CNV and variant calls.
    """
    paired = vcfutils.get_paired_bams([dd.get_align_bam(d) for d in items], items)
    if paired and paired.normal_bam:
        return _get_variants(paired.tumor_data) and _get_cnvs(paired.tumor_data)

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

def estimate(items, batch, config):
    """Estimate heterogeneity for a pair of tumor/normal samples. Run in parallel.

    XXX In progress, currently uses THetA but not yet turned on
    """
    paired = vcfutils.get_paired_bams([dd.get_align_bam(d) for d in items], items)
    cnvs = _get_cnvs(paired.tumor_data)
    new_cnvs = theta.run(cnvs[0], paired)
    print(new_cnvs)

    out = []
    for data in items:
        if batch == _get_batches(data)[0]:
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