"""High level parallel SNP and indel calling using multiple variant callers.
"""
import os
import collections
import copy
import pprint

import six
import toolz as tz

from bcbio import bam, utils
from bcbio.cwl import cwlutils
from bcbio.distributed.split import (grouped_parallel_split_combine, parallel_split_combine)
from bcbio.distributed import multi as dmulti
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import region as pregion
from bcbio.pipeline import shared as pshared
from bcbio.variation import (gatk, gatkfilter, germline, multi,
                             ploidy, vcfutils, vfilter)

# ## Variant filtration -- shared functionality

def variant_filtration(call_file, ref_file, vrn_files, data, items):
    """Filter variant calls using Variant Quality Score Recalibration.

    Newer GATK with Haplotype calling has combined SNP/indel filtering.
    """
    caller = data["config"]["algorithm"].get("variantcaller")
    if "gvcf" not in dd.get_tools_on(data):
        call_file = ploidy.filter_vcf_by_sex(call_file, items)
    if caller in ["freebayes"]:
        return vfilter.freebayes(call_file, ref_file, vrn_files, data)
    elif caller in ["platypus"]:
        return vfilter.platypus(call_file, data)
    elif caller in ["samtools"]:
        return vfilter.samtools(call_file, data)
    elif caller in ["gatk", "gatk-haplotype", "haplotyper"]:
        if dd.get_analysis(data).lower().find("rna-seq") >= 0:
            from bcbio.rnaseq import variation as rnaseq_variation
            return rnaseq_variation.gatk_filter_rnaseq(call_file, data)
        else:
            return gatkfilter.run(call_file, ref_file, vrn_files, data)
    # no additional filtration for callers that filter as part of call process
    else:
        return call_file

# ## High level functionality to run genotyping in parallel

def get_variantcaller(data, key="variantcaller", default=None, require_bam=True):
    if not require_bam or data.get("align_bam"):
        return tz.get_in(["config", "algorithm", key], data, default)

def combine_multiple_callers(samples):
    """Collapse together variant calls from multiple approaches into single data item with `variants`.
    """
    by_bam = collections.OrderedDict()
    for data in (x[0] for x in samples):
        work_bam = tz.get_in(("combine", "work_bam", "out"), data, data.get("align_bam"))
        # For pre-computed VCF inputs, we don't have BAM files
        if not work_bam:
            work_bam = dd.get_sample_name(data)
        jointcaller = tz.get_in(("config", "algorithm", "jointcaller"), data)
        variantcaller = get_variantcaller(data)
        key = (multi.get_batch_for_key(data), work_bam)
        if key not in by_bam:
            by_bam[key] = []
        by_bam[key].append((variantcaller, jointcaller, data))
    out = []
    for callgroup in by_bam.values():
        ready_calls = []
        for variantcaller, jointcaller, data in callgroup:
            if variantcaller:
                cur = data.get("vrn_file_plus", {})
                cur.update({"variantcaller": variantcaller,
                            "vrn_file": data.get("vrn_file_orig") if jointcaller else data.get("vrn_file"),
                            "vrn_file_batch": data.get("vrn_file_batch") if not jointcaller else None,
                            "vrn_stats": data.get("vrn_stats"),
                            "validate": data.get("validate") if not jointcaller else None})
                if jointcaller:
                    cur["population"] = False
                ready_calls.append(cur)
            if jointcaller:
                cur = {"variantcaller": jointcaller,
                       "vrn_file": data.get("vrn_file"),
                       "vrn_file_batch": data.get("vrn_file_batch"),
                       "validate": data.get("validate"),
                       "do_upload": False}
                if not variantcaller:
                    cur["population"] = {"vcf": data.get("vrn_file")}
                ready_calls.append(cur)
            if not jointcaller and not variantcaller:
                ready_calls.append({"variantcaller": "precalled",
                                    "vrn_file": data.get("vrn_file"),
                                    "validate": data.get("validate"),
                                    "do_upload": False})
        final = callgroup[0][-1]
        def orig_variantcaller_order(x):
            try:
                return final["config"]["algorithm"]["orig_variantcaller"].index(x["variantcaller"])
            except ValueError:
                return final["config"]["algorithm"]["orig_jointcaller"].index(x["variantcaller"])
        if len(ready_calls) > 1 and "orig_variantcaller" in final["config"]["algorithm"]:
            final["variants"] = sorted(ready_calls, key=orig_variantcaller_order)
            final["config"]["algorithm"]["variantcaller"] = final["config"]["algorithm"].pop("orig_variantcaller")
            if "orig_jointcaller" in final["config"]["algorithm"]:
                final["config"]["algorithm"]["jointcaller"] = final["config"]["algorithm"].pop("orig_jointcaller")
        else:
            final["variants"] = ready_calls
        final.pop("vrn_file_batch", None)
        final.pop("vrn_file_orig", None)
        final.pop("vrn_file_plus", None)
        final.pop("vrn_stats", None)
        out.append([final])
    return out

def _split_by_ready_regions(ext, file_key, dir_ext_fn):
    """Organize splits based on regions generated by parallel_prep_region.

    Sort splits so largest regions analyzed first, avoiding potentially lagging runs
    at end.
    """
    def _sort_by_size(region_w_bams):
        region, _ = region_w_bams
        _, start, end = region
        return end - start
    def _assign_bams_to_regions(data):
        """Ensure BAMs aligned with input regions, either global or individual.
        """
        for i, region in enumerate(data["region"]):
            work_bams = []
            for xs in data["region_bams"]:
                if len(xs) == 1:
                    work_bams.append(xs[0])
                else:
                    work_bams.append(xs[i])
            for work_bam in work_bams:
                assert os.path.exists(work_bam), work_bam
            yield region, work_bams
    def _do_work(data):
        if "region" in data:
            name = data["group"][0] if "group" in data else data["description"]
            out_dir = os.path.join(data["dirs"]["work"], dir_ext_fn(data))
            out_file = os.path.join(out_dir, "%s%s" % (name, ext))
            assert isinstance(data["region"], (list, tuple))
            out_parts = []
            for r, work_bams in sorted(_assign_bams_to_regions(data), key=_sort_by_size, reverse=True):
                out_region_dir = os.path.join(out_dir, r[0])
                out_region_file = os.path.join(out_region_dir,
                                               "%s-%s%s" % (name, pregion.to_safestr(r), ext))
                out_parts.append((r, work_bams, out_region_file))
            return out_file, out_parts
        else:
            return None, []
    return _do_work

def _collapse_by_bam_variantcaller(samples):
    """Collapse regions to a single representative by BAM input, variant caller and batch.
    """
    by_bam = collections.OrderedDict()
    for data in (x[0] for x in samples):
        work_bam = utils.get_in(data, ("combine", "work_bam", "out"), data.get("align_bam"))
        variantcaller = get_variantcaller(data)
        if isinstance(work_bam, list):
            work_bam = tuple(work_bam)
        key = (multi.get_batch_for_key(data), work_bam, variantcaller)
        try:
            by_bam[key].append(data)
        except KeyError:
            by_bam[key] = [data]
    out = []
    for grouped_data in by_bam.values():
        cur = grouped_data[0]
        cur.pop("region", None)
        region_bams = cur.pop("region_bams", None)
        if region_bams and len(region_bams[0]) > 1:
            cur.pop("work_bam", None)
        out.append([cur])
    return out

def _dup_samples_by_variantcaller(samples, require_bam=True):
    """Prepare samples by variant callers, duplicating any with multiple callers.
    """
    samples = [utils.to_single_data(x) for x in samples]
    samples = germline.split_somatic(samples)
    to_process = []
    extras = []
    for data in samples:
        added = False
        for i, add in enumerate(handle_multiple_callers(data, "variantcaller", require_bam=require_bam)):
            added = True
            add = dd.set_variantcaller_order(add, i)
            to_process.append([add])
        if not added:
            data = _handle_precalled(data)
            data = dd.set_variantcaller_order(data, 0)
            extras.append([data])
    return to_process, extras

def parallel_variantcall_region(samples, run_parallel):
    """Perform variant calling and post-analysis on samples by region.
    """
    to_process, extras = _dup_samples_by_variantcaller(samples)
    split_fn = _split_by_ready_regions(".vcf.gz", "work_bam", get_variantcaller)
    samples = _collapse_by_bam_variantcaller(
        grouped_parallel_split_combine(to_process, split_fn,
                                       multi.group_batches, run_parallel,
                                       "variantcall_sample", "concat_variant_files",
                                       "vrn_file", ["region", "sam_ref", "config"]))
    return extras + samples


def vc_output_record(samples):
    """Prepare output record from variant calling to feed into downstream analysis.

    Prep work handles reformatting so we return generated dictionaries.

    For any shared keys that are calculated only once for a batch, like variant calls
    for the batch, we assign to every sample.
    """
    shared_keys = [["vrn_file"], ["validate", "summary"],
                   ["validate", "tp"], ["validate", "fp"], ["validate", "fn"]]
    raw = cwlutils.samples_to_records([utils.to_single_data(x) for x in samples])
    shared = {}
    for key in shared_keys:
        cur = list(set([x for x in [tz.get_in(key, d) for d in raw] if x]))
        if len(cur) > 0:
            assert len(cur) == 1, (key, cur)
            shared[tuple(key)] = cur[0]
        else:
            shared[tuple(key)] = None
    out = []
    for d in raw:
        for key, val in shared.items():
            d = tz.update_in(d, key, lambda x: val)
        out.append([d])
    return out

def is_joint(data):
    return "gvcf" in dd.get_tools_on(data) or dd.get_jointcaller(data)

def batch_for_variantcall(samples):
    """Prepare a set of samples for parallel variant calling.

    CWL input target that groups samples into batches and variant callers
    for parallel processing.

    If doing joint calling, with `tools_on: [gvcf]`, split the sample into
    individuals instead of combining into a batch.
    """
    sample_order = [dd.get_sample_name(utils.to_single_data(x)) for x in samples]
    to_process, extras = _dup_samples_by_variantcaller(samples, require_bam=False)
    batch_groups = collections.defaultdict(list)
    to_process = [utils.to_single_data(x) for x in to_process]
    for data in cwlutils.samples_to_records(to_process):
        vc = get_variantcaller(data, require_bam=False)
        batches = dd.get_batches(data) or dd.get_sample_name(data)
        if not isinstance(batches, (list, tuple)):
            batches = [batches]
        for b in batches:
            batch_groups[(b, vc)].append(utils.deepish_copy(data))
    batches = []
    for cur_group in batch_groups.values():
        joint_calling = any([is_joint(d) for d in cur_group])
        if joint_calling:
            for d in cur_group:
                batches.append([d])
        else:
            batches.append(cur_group)
    def by_original_order(xs):
        return (min([sample_order.index(dd.get_sample_name(x)) for x in xs]),
                min([dd.get_variantcaller_order(x) for x in xs]))
    return sorted(batches + extras, key=by_original_order)

def _handle_precalled(data):
    """Copy in external pre-called variants fed into analysis.

    Symlinks for non-CWL runs where we want to ensure VCF present
    in a local directory.
    """
    if data.get("vrn_file") and not cwlutils.is_cwl_run(data):
        vrn_file = data["vrn_file"]
        if isinstance(vrn_file, (list, tuple)):
            assert len(vrn_file) == 1
            vrn_file = vrn_file[0]
        precalled_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "precalled"))
        ext = utils.splitext_plus(vrn_file)[-1]
        orig_file = os.path.abspath(vrn_file)
        our_vrn_file = os.path.join(precalled_dir, "%s-precalled%s" % (dd.get_sample_name(data), ext))
        utils.copy_plus(orig_file, our_vrn_file)
        data["vrn_file"] = our_vrn_file
    return data

def handle_multiple_callers(data, key, default=None, require_bam=True):
    """Split samples that potentially require multiple variant calling approaches.
    """
    callers = get_variantcaller(data, key, default, require_bam=require_bam)
    if isinstance(callers, six.string_types):
        return [data]
    elif not callers:
        return []
    else:
        out = []
        for caller in callers:
            base = copy.deepcopy(data)
            if not base["config"]["algorithm"].get("orig_%s" % key):
                base["config"]["algorithm"]["orig_%s" % key] = \
                  base["config"]["algorithm"][key]
            base["config"]["algorithm"][key] = caller
            # if splitting by variant caller, also split by jointcaller
            if key == "variantcaller":
                jcallers = get_variantcaller(data, "jointcaller", [])
                if isinstance(jcallers, six.string_types):
                    jcallers = [jcallers]
                if jcallers:
                    base["config"]["algorithm"]["orig_jointcaller"] = jcallers
                    jcallers = [x for x in jcallers if x.startswith(caller)]
                    if jcallers:
                        base["config"]["algorithm"]["jointcaller"] = jcallers[0]
                    else:
                        base["config"]["algorithm"]["jointcaller"] = False
            out.append(base)
        return out


# Avoid use of multicore GATK HaplotypeCallerSpark due to bugs in gVCF output
# during joint calling Intermittent failures where we don't produce correct headers
# https://github.com/broadinstitute/gatk/issues/4821
# "gatk-haplotype"
SUPPORT_MULTICORE = ["strelka2", "haplotyper", "tnhaplotyper", "tnscope",
                     "deepvariant", "pisces", "octopus", "smcounter2"]

def get_variantcallers():
    from bcbio.variation import (freebayes, cortex, samtools, varscan, mutect, mutect2, octopus,
                                 pisces, platypus, scalpel, sentieon, strelka2, vardict, qsnp,
                                 deepvariant, smcounter2)
    return {"gatk": gatk.unified_genotyper,
            "gatk-haplotype": gatk.haplotype_caller,
            "mutect2": mutect2.mutect2_caller,
            "freebayes": freebayes.run_freebayes,
            "cortex": cortex.run_cortex,
            "deepvariant": deepvariant.run,
            "samtools": samtools.run_samtools,
            "varscan": varscan.run_varscan,
            "mutect": mutect.mutect_caller,
            "octopus": octopus.run,
            "pisces": pisces.run,
            "platypus": platypus.run,
            "scalpel": scalpel.run_scalpel,
            "smcounter2": smcounter2.run,
            "strelka2": strelka2.run,
            "vardict": vardict.run_vardict,
            "vardict-java": vardict.run_vardict,
            "vardict-perl": vardict.run_vardict,
            "haplotyper": sentieon.run_haplotyper,
            "tnhaplotyper": sentieon.run_tnhaplotyper,
            "tnscope": sentieon.run_tnscope,
            "qsnp": qsnp.run_qsnp}

def variantcall_sample(data, region=None, align_bams=None, out_file=None):
    """Parallel entry point for doing genotyping of a region of a sample.
    """
    if out_file is None or not os.path.exists(out_file) or not os.path.lexists(out_file):
        utils.safe_makedir(os.path.dirname(out_file))
        ref_file = dd.get_ref_file(data)
        config = data["config"]
        caller_fns = get_variantcallers()
        caller_fn = caller_fns[config["algorithm"].get("variantcaller")]
        if len(align_bams) == 1:
            items = [data]
        else:
            items = multi.get_orig_items(data)
            assert len(items) == len(align_bams)
        assoc_files = tz.get_in(("genome_resources", "variation"), data, {})
        if not assoc_files: assoc_files = {}
        for bam_file in align_bams:
            bam.index(bam_file, data["config"], check_timestamp=False)
        out_file = caller_fn(align_bams, items, ref_file, assoc_files, region, out_file)
    if region:
        data["region"] = region
    data["vrn_file"] = out_file
    return [data]

def concat_batch_variantcalls(items, region_block=True, skip_jointcheck=False):
    """CWL entry point: combine variant calls from regions into single VCF.
    """
    items = [utils.to_single_data(x) for x in items]
    batch_name = _get_batch_name(items, skip_jointcheck)
    variantcaller = _get_batch_variantcaller(items)
    # Pre-called input variant files
    if not variantcaller and all(d.get("vrn_file") for d in items):
        return {"vrn_file": items[0]["vrn_file"]}
    out_file = os.path.join(dd.get_work_dir(items[0]), variantcaller, "%s.vcf.gz" % (batch_name))
    utils.safe_makedir(os.path.dirname(out_file))
    if region_block:
        regions = [_region_to_coords(rs[0]) for rs in items[0]["region_block"]]
    else:
        regions = [_region_to_coords(r) for r in items[0]["region"]]
    vrn_file_regions = items[0]["vrn_file_region"]
    out_file = vcfutils.concat_variant_files(vrn_file_regions, out_file, regions,
                                             dd.get_ref_file(items[0]), items[0]["config"])
    return {"vrn_file": out_file}

def _region_to_coords(region):
    """Split GATK region specification (chr1:1-10) into a tuple of chrom, start, end
    """
    chrom, coords = region.split(":")
    start, end = coords.split("-")
    return (chrom, int(start), int(end))

def _get_batch_name(items, skip_jointcheck=False):
    """Retrieve the shared batch name for a group of items.
    """
    batch_names = collections.defaultdict(int)
    has_joint = any([is_joint(d) for d in items])
    for data in items:
        if has_joint and not skip_jointcheck:
            batches = dd.get_sample_name(data)
        else:
            batches = dd.get_batches(data) or dd.get_sample_name(data)
        if not isinstance(batches, (list, tuple)):
            batches = [batches]
        for b in batches:
            batch_names[b] += 1
    return sorted(batch_names.items(), key=lambda x: x[-1], reverse=True)[0][0]

def _get_batch_variantcaller(items):
    variantcaller = [vc for vc in list(set([get_variantcaller(x, require_bam=False) for x in items])) if vc]
    if not variantcaller:
        return None
    assert len(variantcaller) == 1, "%s\n%s" % (variantcaller, pprint.pformat(items))
    return variantcaller[0]

def variantcall_batch_region(items):
    """CWL entry point: variant call a batch of samples in a block of regions.
    """
    items = [utils.to_single_data(x) for x in items]
    align_bams = [dd.get_align_bam(x) for x in items]
    variantcaller = _get_batch_variantcaller(items)
    region_blocks = list(set([tuple(x.get("region_block")) for x in items if "region_block" in x]))
    assert len(region_blocks) == 1, region_blocks
    region_block = region_blocks[0]
    # Pre-called input variant files
    if not variantcaller and all(d.get("vrn_file") for d in items):
        return {"vrn_file_region": None, "region_block": region_block}
    caller_fn = get_variantcallers()[variantcaller]
    assoc_files = tz.get_in(("genome_resources", "variation"), items[0], {})
    region = _region_to_coords(region_block[0])
    chrom, start, end = region
    region_str = "_".join(str(x) for x in region)
    batch_name = _get_batch_name(items)
    out_file = os.path.join(dd.get_work_dir(items[0]), variantcaller, chrom,
                            "%s-%s-block.vcf.gz" % (batch_name, region_str))
    utils.safe_makedir(os.path.dirname(out_file))
    with pshared.bedtools_tmpdir(items[0]):
        if variantcaller in SUPPORT_MULTICORE:
            call_file = caller_fn(align_bams, items, dd.get_ref_file(items[0]), assoc_files,
                                  [_region_to_coords(r) for r in region_block], out_file)
        else:
            call_file = _run_variantcall_batch_multicore(items, region_block, out_file)
    return {"vrn_file_region": call_file, "region_block": region_block}

def _run_variantcall_batch_multicore(items, regions, final_file):
    """Run variant calling on a batch of items using multiple cores.
    """
    batch_name = _get_batch_name(items)
    variantcaller = _get_batch_variantcaller(items)
    work_bams = [dd.get_work_bam(d) or dd.get_align_bam(d) for d in items]
    def split_fn(data):
        out = []
        for region in regions:
            region = _region_to_coords(region)
            chrom, start, end = region
            region_str = "_".join(str(x) for x in region)
            out_file = os.path.join(dd.get_work_dir(items[0]), variantcaller, chrom,
                                    "%s-%s.vcf.gz" % (batch_name, region_str))
            out.append((region, work_bams, out_file))
        return final_file, out
    parallel = {"type": "local", "num_jobs": dd.get_num_cores(items[0]), "cores_per_job": 1}
    run_parallel = dmulti.runner(parallel, items[0]["config"])
    to_run = copy.deepcopy(items[0])
    to_run["sam_ref"] = dd.get_ref_file(to_run)
    to_run["group_orig"] = items
    parallel_split_combine([[to_run]], split_fn, run_parallel,
                           "variantcall_sample", "concat_variant_files",
                           "vrn_file", ["region", "sam_ref", "config"])
    return final_file
