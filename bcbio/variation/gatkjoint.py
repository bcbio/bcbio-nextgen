"""Perform joint genotyping using GATK HaplotypeCaller with gVCF inputs

Handles merging of large batch sizes using CombineGVCFs and
joint variant calling with GenotypeGVCFs.
"""
import math
import toolz as tz

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.variation import bamprep, vcfutils

def run_region(data, region, vrn_files, out_file):
    """Perform variant calling on gVCF inputs in a specific genomic region.
    """
    vrn_files = _batch_gvcfs(data, region, vrn_files, dd.get_ref_file(data), out_file)
    return _run_genotype_gvcfs(data, region, vrn_files, dd.get_ref_file(data), out_file)

# ## gVCF joint genotype calling

def _run_genotype_gvcfs(data, region, vrn_files, ref_file, out_file):
    if not utils.file_exists(out_file):
        broad_runner = broad.runner_from_config(data["config"])
        with file_transaction(data, out_file) as tx_out_file:
            assoc_files = tz.get_in(("genome_resources", "variation"), data, {})
            if not assoc_files: assoc_files = {}
            params = ["-T", "GenotypeGVCFs",
                      "-R", ref_file, "-o", tx_out_file,
                      "-L", bamprep.region_to_gatk(region)]
            for vrn_file in vrn_files:
                params += ["--variant", vrn_file]
            if assoc_files.get("dbsnp"):
                params += ["--dbsnp", assoc_files["dbsnp"]]
            broad_runner.new_resources("gatk-haplotype")
            cores = dd.get_cores(data)
            if cores > 1:
                params += ["-nt", str(cores)]
                memscale = {"magnitude": 0.9 * cores, "direction": "increase"}
            else:
                memscale = None
            broad_runner.run_gatk(params, memscale=memscale)
    return vcfutils.bgzip_and_index(out_file, data["config"])

# ## gVCF batching

def _batch_gvcfs(data, region, vrn_files, ref_file, out_file=None):
    """Perform batching of gVCF files if above recommended input count.
    """
    if out_file is None:
        out_file = vrn_files[0]
    # group to get below the maximum batch size, using 200 as the baseline
    max_batch = int(dd.get_joint_group_size(data))
    if len(vrn_files) > max_batch:
        out = []
        num_batches = int(math.ceil(float(len(vrn_files)) / max_batch))
        for i, batch_vrn_files in enumerate(tz.partition_all(num_batches, vrn_files)):
            base, ext = utils.splitext_plus(out_file)
            batch_out_file = "%s-b%s%s" % (base, i, ext)
            out.append(run_combine_gvcfs(batch_vrn_files, region, ref_file, batch_out_file, data))
        return _batch_gvcfs(data, region, out, ref_file)
    else:
        return vrn_files

def run_combine_gvcfs(vrn_files, region, ref_file, out_file, data):
    if not utils.file_exists(out_file):
        broad_runner = broad.runner_from_config(data["config"])
        with file_transaction(data, out_file) as tx_out_file:
            params = ["-T", "CombineGVCFs", "-R", ref_file, "-o", tx_out_file]
            if region:
                params += ["-L", bamprep.region_to_gatk(region)]
            for vrn_file in vrn_files:
                params += ["--variant", vrn_file]
            cores = dd.get_cores(data)
            memscale = {"magnitude": 0.9 * cores, "direction": "increase"} if cores > 1 else None
            broad_runner.new_resources("gatk-haplotype")
            broad_runner.run_gatk(params, memscale=memscale)
    return vcfutils.bgzip_and_index(out_file, data["config"])
