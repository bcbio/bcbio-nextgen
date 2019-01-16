"""Detect structural variation in genomes using high-throughput sequencing data.
"""
import collections
import copy
import operator
import os

import toolz as tz

from bcbio import utils
from bcbio.cwl import cwlutils
from bcbio.pipeline import datadict as dd
from bcbio.structural import (battenberg, cn_mops, cnvkit, delly, gatkcnv, gridss,
                              lumpy, manta, metasv, prioritize, purecn, purple, plot,
                              seq2c, titancna, validate, wham)
from bcbio.variation import validate as vcvalidate
from bcbio.variation import vcfutils

import six
from functools import reduce


# Stratify callers by stage -- see `run` documentation below for definitions
_CALLERS = {
  "initial": {"cnvkit": cnvkit.run},
  "standard": {"cn.mops": cn_mops.run, "manta": manta.run, "cnvkit": cnvkit.run,
               "delly": delly.run, "lumpy": lumpy.run, "wham": wham.run,
               "battenberg": battenberg.run, "seq2c": seq2c.run, "gridss": gridss.run,
               "titancna": titancna.run, "purecn": purecn.run, "purple": purple.run,
               "gatk-cnv": gatkcnv.run},
  "ensemble": {"metasv": metasv.run,
               "prioritize": prioritize.run}}
_NEEDS_BACKGROUND = set(["cn.mops"])
_GLOBAL_BATCHING = set(["seq2c"])
# CNV callers that have background references
_CNV_REFERENCE = set(["seq2c", "cnvkit", "gatk-cnv"])

def _get_callers(items, stage, special_cases=False):
    """Retrieve available callers for the provided stage.

    Handles special cases like CNVkit that can be in initial or standard
    depending on if fed into Lumpy analysis.
    """
    callers = utils.deepish_copy(_CALLERS[stage])
    if special_cases and "cnvkit" in callers:
        has_lumpy = any("lumpy" in get_svcallers(d) or "lumpy" in d["config"]["algorithm"].get("svcaller_orig", [])
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
    elif isinstance(svs, six.string_types):
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
            # clean SV callers present in multiple rounds and not this caller
            final_svs = []
            for sv in data.get("sv", []):
                if (stage == "ensemble" or sv["variantcaller"] == svcaller or sv["variantcaller"] not in svs
                      or svcaller not in _get_callers([data], stage, special_cases=True)):
                    final_svs.append(sv)
            base["sv"] = final_svs
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
        if len(batches) > 1:
            lead_batches[(dd.get_sample_name(final), dd.get_phenotype(final) == "germline")] = batches[0]
        for batch in batches:
            try:
                by_batch[batch].append(final)
            except KeyError:
                by_batch[batch] = [final]
    out = []
    for batch, items in by_batch.items():
        if any("svplots" in dd.get_tools_on(d) for d in items):
            items = plot.by_regions(items)
        for data in items:
            if lead_batches.get((dd.get_sample_name(data), dd.get_phenotype(data) == "germline")) in [batch, None]:
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
    samples = cwlutils.assign_complex_to_samples(samples)
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
                if stage in ["ensemble"]:  # no batching for ensemble methods
                    if isinstance(batch, six.string_types) and batch != dd.get_sample_name(x):
                        batch += "_%s" % dd.get_sample_name(x)
                    else:
                        batch = dd.get_sample_name(x)
                    if dd.get_phenotype(x) == "germline":
                        batch += "_germline"
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
    items = [utils.to_single_data(x) for x in items]
    items = cwlutils.unpack_tarballs(items, items[0])
    svcaller = items[0]["config"]["algorithm"].get("svcaller")
    caller_fn = _get_callers(items, stage, special_cases=True).get(svcaller)
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
    # Avoid nesting of callers for CWL runs for easier extraction
    if cwlutils.is_cwl_run(items[0]):
        out_cwl = []
        for data in [utils.to_single_data(x) for x in out]:
            # Run validation directly from CWL runs since we're single stage
            data = validate.evaluate(data)
            data["svvalidate"] = {"summary": tz.get_in(["sv-validate", "csv"], data)}
            svs = data.get("sv")
            if svs:
                assert len(svs) == 1, svs
                data["sv"] = svs[0]
            else:
                data["sv"] = {}
            data = _add_supplemental(data)
            out_cwl.append([data])
        return out_cwl
    return out

def _add_supplemental(data):
    """Add additional supplemental files to CWL sv output, give useful names.
    """
    if "supplemental" not in data["sv"]:
        data["sv"]["supplemental"] = []
    if data["sv"].get("variantcaller"):
        cur_name = _useful_basename(data)
        for k in ["cns", "vrn_bed"]:
            if data["sv"].get(k) and os.path.exists(data["sv"][k]):
                dname, orig = os.path.split(data["sv"][k])
                orig_base, orig_ext = utils.splitext_plus(orig)
                orig_base = _clean_name(orig_base, data)
                if orig_base:
                    fname = "%s-%s%s" % (cur_name, orig_base, orig_ext)
                else:
                    fname = "%s%s" % (cur_name, orig_ext)
                sup_out_file = os.path.join(dname, fname)
                utils.symlink_plus(data["sv"][k], sup_out_file)
                data["sv"]["supplemental"].append(sup_out_file)
    return data

def _clean_name(fname, data):
    """Remove standard prefixes from a filename before renaming with useful names.
    """
    for to_remove in dd.get_batches(data) + [dd.get_sample_name(data), data["sv"]["variantcaller"]]:
        for ext in ("-", "_"):
            if fname.startswith("%s%s" % (to_remove, ext)):
                fname = fname[len(to_remove) + len(ext):]
        if fname.startswith(to_remove):
            fname = fname[len(to_remove):]
    return fname

def _useful_basename(data):
    """Provide a useful file basename for outputs, referencing batch/sample and caller.
    """
    names = dd.get_batches(data)
    if not names:
        names = [dd.get_sample_name(data)]
    batch_name = names[0]
    return "%s-%s" % (batch_name, data["sv"]["variantcaller"])

def _group_by_sample(items):
    """Group a set of items by sample names + multiple callers for prioritization
    """
    by_sample = collections.defaultdict(list)
    for d in items:
        by_sample[dd.get_sample_name(d)].append(d)
    out = []
    for sample_group in by_sample.values():
        cur = utils.deepish_copy(sample_group[0])
        svs = []
        for d in sample_group:
            svs.append(d["sv"])
        cur["sv"] = svs
        out.append(cur)
    return out

def summarize_sv(items):
    """CWL target: summarize structural variants for multiple samples.

    XXX Need to support non-VCF output as tabix indexed output
    """
    items = [utils.to_single_data(x) for x in vcvalidate.summarize_grading(items, "svvalidate")]
    out = {"sv": {"calls": [],
                  "supplemental": [],
                  "prioritize": {"tsv": [],
                                 "raw": []}},
           "svvalidate": vcvalidate.combine_validations(items, "svvalidate")}
    added = set([])
    # Standard callers
    for data in items:
        if data.get("sv"):
            if data["sv"].get("vrn_file"):
                ext = utils.splitext_plus(data["sv"]["vrn_file"])[-1]
                cur_name = _useful_basename(data)
                if cur_name not in added and ext.startswith(".vcf"):
                    added.add(cur_name)
                    out_file = os.path.join(utils.safe_makedir(os.path.join(dd.get_work_dir(data),
                                                                            "sv", "calls")),
                                            "%s%s" % (cur_name, ext))
                    utils.copy_plus(data["sv"]["vrn_file"], out_file)
                    out_file = vcfutils.bgzip_and_index(out_file, data["config"])
                    out["sv"]["calls"].append(out_file)
            if data["sv"].get("supplemental"):
                out["sv"]["supplemental"].extend([x for x in data["sv"]["supplemental"] if x])
    # prioritization
    for pdata in _group_by_sample(items):
        prioritysv = [x for x in prioritize.run([utils.deepish_copy(pdata)])[0].get("sv", [])
                      if x["variantcaller"] == "sv-prioritize"]
        if prioritysv:
            out["sv"]["prioritize"]["tsv"].append(prioritysv[0]["vrn_file"])
            out["sv"]["prioritize"]["raw"].extend(prioritysv[0]["raw_files"].values())
    return [out]

# ## configuration

def standardize_cnv_reference(data):
    """Standardize cnv_reference background to support multiple callers.
    """
    out = tz.get_in(["config", "algorithm", "background", "cnv_reference"], data, {})
    cur_callers = set(data["config"]["algorithm"].get("svcaller")) & _CNV_REFERENCE
    if isinstance(out, six.string_types):
        if not len(cur_callers) == 1:
            raise ValueError("Multiple CNV callers and single background reference for %s: %s" %
                                data["description"], list(cur_callers))
        else:
            out = {cur_callers.pop(): out}
    return out

def supports_cnv_reference(c):
    return c in _CNV_REFERENCE

def parallel_multiplier(items):
    """Use more resources (up to available limits) if we have multiple QC samples/svcallers.
    """
    machines = []
    for data in (xs[0] for xs in items):
        machines.append(max(1, len(get_svcallers(data)), len(dd.get_algorithm_qc(data))))
    return sum(machines)
