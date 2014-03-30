"""High level parallel SNP and indel calling using multiple variant callers.
"""
import os
import collections
import copy

from bcbio import utils
from bcbio.distributed.split import grouped_parallel_split_combine
from bcbio.pipeline.shared import process_bam_by_chromosome
from bcbio.variation import gatk, gatkfilter, multi, phasing, ploidy, vfilter

# ## Variant filtration -- shared functionality

def variant_filtration(call_file, ref_file, vrn_files, data):
    """Filter variant calls using Variant Quality Score Recalibration.

    Newer GATK with Haplotype calling has combined SNP/indel filtering.
    """
    caller = data["config"]["algorithm"].get("variantcaller")
    call_file = ploidy.filter_vcf_by_sex(call_file, data)
    if caller in ["freebayes"]:
        return vfilter.freebayes(call_file, ref_file, vrn_files, data)
    elif caller in ["gatk", "gatk-haplotype"]:
        return gatkfilter.run(call_file, ref_file, vrn_files, data)
    # no additional filtration for callers that filter as part of call process
    else:
        return call_file

# ## High level functionality to run genotyping in parallel

def get_variantcaller(data):
    return data["config"]["algorithm"].get("variantcaller", "gatk")

def combine_multiple_callers(data):
    """Collapse together variant calls from multiple approaches into variants
    """
    by_bam = collections.defaultdict(list)
    for x in data:
        by_bam[x[0]["work_bam"]].append(x[0])
    out = []
    for grouped_calls in by_bam.itervalues():
        ready_calls = [{"variantcaller": get_variantcaller(x),
                        "vrn_file": x.get("vrn_file"),
                        "validate": x.get("validate")}
                       for x in grouped_calls]
        final = grouped_calls[0]
        def orig_variantcaller_order(x):
            return final["config"]["algorithm"]["orig_variantcaller"].index(x["variantcaller"])
        if len(ready_calls) > 1 and "orig_variantcaller" in final["config"]["algorithm"]:
            final["variants"] = sorted(ready_calls, key=orig_variantcaller_order)
        else:
            final["variants"] = ready_calls
        out.append([final])
    return out

def handle_multiple_variantcallers(data):
    """Split samples that potentially require multiple variant calling approaches.
    """
    assert len(data) == 1
    callers = get_variantcaller(data[0])
    if isinstance(callers, basestring):
        return [data]
    elif not callers:
        return []
    else:
        out = []
        for caller in callers:
            base = copy.deepcopy(data[0])
            base["config"]["algorithm"]["orig_variantcaller"] = \
              base["config"]["algorithm"]["variantcaller"]
            base["config"]["algorithm"]["variantcaller"] = caller
            out.append([base])
        return out

def parallel_variantcall(sample_info, parallel_fn):
    """Provide sample genotyping, running in parallel over individual chromosomes.
    """
    to_process = []
    finished = []
    for x in sample_info:
        if get_variantcaller(x[0]):
            to_process.extend(handle_multiple_variantcallers(x))
        else:
            finished.append(x)
    if len(to_process) > 0:
        split_fn = process_bam_by_chromosome("-variants.vcf.gz", "work_bam",
                                             dir_ext_fn=get_variantcaller)
        processed = grouped_parallel_split_combine(
            to_process, split_fn, multi.group_batches, parallel_fn,
            "variantcall_sample", "split_variants_by_sample", "combine_variant_files",
            "vrn_file", ["sam_ref", "config"])
        finished.extend(processed)
    return finished

def get_variantcallers():
    from bcbio.variation import freebayes, cortex, samtools, varscan, mutect
    return {"gatk": gatk.unified_genotyper,
            "gatk-haplotype": gatk.haplotype_caller,
            "freebayes": freebayes.run_freebayes,
            "cortex": cortex.run_cortex,
            "samtools": samtools.run_samtools,
            "varscan": varscan.run_varscan,
            "mutect": mutect.mutect_caller}

def variantcall_sample(data, region=None, out_file=None):
    """Parallel entry point for doing genotyping of a region of a sample.
    """
    if out_file is None or not os.path.exists(out_file) or not os.path.lexists(out_file):
        utils.safe_makedir(os.path.dirname(out_file))
        sam_ref = data["sam_ref"]
        config = data["config"]
        caller_fns = get_variantcallers()
        caller_fn = caller_fns[config["algorithm"].get("variantcaller", "gatk")]
        if isinstance(data["work_bam"], basestring):
            align_bams = [data["work_bam"]]
            items = [data]
        else:
            align_bams = data["work_bam"]
            items = data["work_items"]
        call_file = "%s-raw%s" % utils.splitext_plus(out_file)
        call_file = caller_fn(align_bams, items, sam_ref,
                              data["genome_resources"]["variation"],
                              region, call_file)
        if data["config"]["algorithm"].get("phasing", False) == "gatk":
            call_file = phasing.read_backed_phasing(call_file, align_bams, sam_ref, region, config)
        utils.symlink_plus(call_file, out_file)
        if "work_items" in data:
            del data["work_items"]
    data["vrn_file"] = out_file
    return [data]
