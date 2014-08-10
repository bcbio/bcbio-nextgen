"""Perform joint genotyping using GATK HaplotypeCaller with gVCF inputs

Handles merging of large batch sizes using CombineGVCFs and
joint variant calling with GenotypeGVCFs.
"""

import toolz as tz

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.variation import bamprep

def run_region(data, region, vrn_files, out_file):
    """Perform variant calling on gVCF inputs in a specific genomic region.
    """
    vrn_files = _batch_gvcfs(data, region, vrn_files, dd.get_ref_file(data), out_file)
    return _run_genotype_gvcfs(data, region, vrn_files, dd.get_ref_file(data), out_file)

# ## gVCF joint genotype calling

def _run_genotype_gvcfs(data, region, vrn_files, ref_file, out_file):
    if not utils.file_exists(out_file):
        broad_runner = broad.runner_from_config(data["config"])
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "GenotypeGVCFs", "-R", ref_file, "-o", tx_out_file,
                      "-L", bamprep.region_to_gatk(region)]
            for vrn_file in vrn_files:
                params += ["--variant", vrn_file]
            broad_runner.new_resources("gatk-haplotype")
            broad_runner.run_gatk(params)
    return out_file

# ## gVCF batching

MAX_BATCH = 200  # Recommended sample count from GATK team where we should use CombineGVCFs

def _batch_gvcfs(data, region, vrn_files, ref_file, out_file=None):
    """Perform batching of gVCF files if above recommended input count.
    """
    if out_file is None:
        out_file = vrn_files[0]
    if len(vrn_files) >= MAX_BATCH:
        out = []
        for i, batch_vrn_files in enumerate(tz.partition_all(MAX_BATCH, vrn_files)):
            base, ext = utils.splitext_plus(out_file)
            batch_out_file = "%s-b%s%s" % (base, i, ext)
            out.append(_run_combine_gvcfs(batch_vrn_files, region, ref_file, batch_out_file, data))
        return _batch_gvcfs(data, region, out, ref_file)
    else:
        return vrn_files

def _run_combine_gvcfs(vrn_files, region, ref_file, out_file, data):
    if not utils.file_exists(out_file):
        broad_runner = broad.runner_from_config(data["config"])
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "CombineGVCFs", "-R", ref_file, "-o", tx_out_file,
                      "-L", bamprep.region_to_gatk(region)]
            for vrn_file in vrn_files:
                params += ["--variant", vrn_file]
            broad_runner.new_resources("gatk-haplotype")
            broad_runner.run_gatk(params)
    return out_file
