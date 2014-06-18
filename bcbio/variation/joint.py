"""Joint variant calling with multiple samples: aka squaring off, or backfilling.

Handles the N+1 problem of variant calling by combining and recalling samples
previously calling individually (or in smaller batches). Recalls at all positions found
variable in any of the input samples within each batch. Takes a general approach supporting
GATK's incremental joint discovery (http://www.broadinstitute.org/gatk/guide/article?id=3893)
and FreeBayes's N+1 approach (https://groups.google.com/d/msg/freebayes/-GK4zI6NsYY/Wpcp8nt_PVMJ)
as implemented in bcbio.variation.recall (https://github.com/chapmanb/bcbio.variation.recall).
"""
import collections
import contextlib
import os

try:
    import pybedtools
except ImportError:
    pybedtools = None
import pysam
import toolz as tz

from bcbio import utils
from bcbio.distributed.split import grouped_parallel_split_combine
from bcbio.pipeline import config_utils, region
from bcbio.provenance import do
from bcbio.variation import bamprep, multi

def _get_callable_regions(data):
    """Retrieve regions to parallelize by from callable regions, variant regions or chromosomes
    """
    callable_files = data.get("callable_regions") or data.get("variant_regions")
    if callable_files:
        assert len(callable_files) == 1
        regions = [(r.chrom, int(r.start), int(r.stop)) for r in pybedtools.BedTool(callable_files[0])]
    else:
        work_bam = list(tz.take(1, filter(lambda x: x.endswith(".bam"), data["work_bams"])))
        if work_bam:
            with contextlib.closing(pysam.Samfile(work_bam[0], "rb")) as pysam_bam:
                regions = [(chrom, 0, length) for (chrom, length) in zip(pysam_bam.references,
                                                                         pysam_bam.lengths)]
        else:
            raise NotImplementedError("No variant regions or BAM files to calculate chromosomes")
    return regions

def _split_by_callable_region(data):
    """Split by callable or variant regions.

    We expect joint calling to be deep in numbers of samples per region, so prefer
    splitting aggressively by regions.
    """
    batch = tz.get_in(("metadata", "batch"), data)
    name = batch if batch else tz.get_in(("rgnames", "sample"), data)
    out_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "joint", name))
    utils.safe_makedir(os.path.join(out_dir, "inprep"))
    parts = []
    for feat in _get_callable_regions(data):
        region_dir = utils.safe_makedir(os.path.join(out_dir, feat[0]))
        region_prep_dir = os.path.join(region_dir, "inprep")
        if not os.path.exists(region_prep_dir):
            os.symlink(os.path.join(os.pardir, "inprep"), region_prep_dir)
        region_outfile = os.path.join(region_dir, "%s-%s.vcf.gz" % (batch, region.to_safestr(feat)))
        parts.append((feat, data["work_bams"], data["vrn_files"], region_outfile))
    out_file = os.path.join(out_dir, "%s-joint.vcf.gz" % name)
    return out_file, parts

def square_off(samples, run_parallel):
    """Perform joint calling at all variants within a batch.
    """
    to_process = []
    extras = []
    for data in [x[0] for x in samples]:
        jointcaller = tz.get_in(("config", "algorithm", "jointcaller"), data)
        batch = tz.get_in(("metadata", "batch"), data)
        if jointcaller and batch:
            to_process.append([data])
        else:
            extras.append([data])
    processed = grouped_parallel_split_combine(to_process, _split_by_callable_region,
                                               multi.group_batches_joint, run_parallel,
                                               "square_batch_region", "concat_variant_files",
                                               "vrn_file", ["region", "sam_ref", "config"])
    return _combine_to_jointcaller(processed) + extras

def _combine_to_jointcaller(processed):
    """Add joint calling information to variants, while collapsing independent regions.
    """
    by_vrn_file = collections.OrderedDict()
    for data in (x[0] for x in processed):
        key = (tz.get_in(("config", "algorithm", "jointcaller"), data), data["vrn_file"])
        if key not in by_vrn_file:
            by_vrn_file[key] = []
        by_vrn_file[key].append(data)
    out = []
    for grouped_data in by_vrn_file.values():
        cur = grouped_data[0]
        out.append([cur])
    return out

def square_batch_region(data, region, bam_files, vrn_files, out_file):
    """Perform squaring of a batch in a supplied region, with input BAMs
    """
    if not utils.file_exists(out_file):
        jointcaller = tz.get_in(("config", "algorithm", "jointcaller"), data)
        if jointcaller == "bcbio-variation-recall":
            _square_batch_bcbio_variation(data, region, bam_files, vrn_files, out_file)
        else:
            raise ValueError("Unexpected joint calling approach: %s" % jointcaller)
    if region:
        data["region"] = region
    data = _fix_orig_vcf_refs(data)
    data["vrn_file"] = out_file
    return [data]

def _fix_orig_vcf_refs(data):
    """Supply references to initial variantcalls if run in addition to batching.
    """
    variantcaller = tz.get_in(("config", "algorithm", "variantcaller"), data)
    if variantcaller:
        data["vrn_file_orig"] = data["vrn_file"]
    for i, sub in enumerate(data["group_orig"]):
        sub_vc = tz.get_in(("config", "algorithm", "variantcaller"), sub)
        if sub_vc:
            sub["vrn_file_orig"] = sub.pop("vrn_file")
            data["group_orig"][i] = sub
    return data

def _square_batch_bcbio_variation(data, region, bam_files, vrn_files, out_file):
    """
    """
    ref_file = tz.get_in(("reference", "fasta", "base"), data)
    cores = tz.get_in(("config", "algorithm", "num_cores"), data, 1)
    resources = config_utils.get_resources("bcbio-variation-recall", data["config"])
    jvm_opts = config_utils.adjust_opts(resources.get("jvm_opts", ["-Xms750m", "-Xmx2g"]),
                                        {"algorithm": {"memory_adjust": {"direction": "increase",
                                                                         "magnitude": cores}}})
    cmd = ["bcbio-variation-recall", "square"] + jvm_opts + \
          ["-c", cores, "-r", bamprep.region_to_gatk(region)] + \
          [out_file, ref_file] + vrn_files + bam_files
    do.run(cmd, "Squaring off in region: %s" % bamprep.region_to_gatk(region))
    return out_file
