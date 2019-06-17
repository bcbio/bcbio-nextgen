"""Perform joint genotyping using GATK HaplotypeCaller with gVCF inputs.

For GATK4, merges into a shared database using GenomicsDBImport. For GATK3
handles merging of large batch sizes using CombineGVCFs.
For both, follows this with joint variant calling using GenotypeGVCFs.
"""
import math
import os
import shutil
import toolz as tz

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.variation import bamprep, ploidy, vcfutils

def run_region(data, region, vrn_files, out_file):
    """Perform variant calling on gVCF inputs in a specific genomic region.
    """
    broad_runner = broad.runner_from_config(data["config"])
    if broad_runner.gatk_type() == "gatk4":
        genomics_db = _run_genomicsdb_import(vrn_files, region, out_file, data)
        return _run_genotype_gvcfs_genomicsdb(genomics_db, region, out_file, data)
    else:
        vrn_files = _batch_gvcfs(data, region, vrn_files, dd.get_ref_file(data), out_file)
        return _run_genotype_gvcfs_gatk3(data, region, vrn_files, dd.get_ref_file(data), out_file)

# ## gVCF joint calling -- GATK4

def _run_genomicsdb_import(vrn_files, region, out_file, data):
    """Create a GenomicsDB reference for all the variation files: GATK4.

    Not yet tested as scale, need to explore --batchSize to reduce memory
    usage if needed.

    Does not support transactional directories yet, since
    GenomicsDB databases cannot be moved to new locations. We try to
    identify half-finished databases and restart:
https://gatkforums.broadinstitute.org/gatk/discussion/10061/using-genomicsdbimport-to-prepare-gvcfs-for-input-to-genotypegvcfs-in-gatk4

    Known issue -- Genomics DB workspace path core dumps on longer paths:
    (std::string::compare(char const*))
    """
    out_dir = "%s_genomicsdb" % utils.splitext_plus(out_file)[0]
    if not os.path.exists(out_dir) or _incomplete_genomicsdb(out_dir):
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        with utils.chdir(os.path.dirname(out_file)):
            with file_transaction(data, out_dir) as tx_out_dir:
                broad_runner = broad.runner_from_config(data["config"])
                cores = dd.get_cores(data)
                params = ["-T", "GenomicsDBImport",
                          "--reader-threads", str(cores),
                          "--genomicsdb-workspace-path", os.path.relpath(out_dir, os.getcwd()),
                          "-L", bamprep.region_to_gatk(region)]
                for vrn_file in vrn_files:
                    vcfutils.bgzip_and_index(vrn_file, data["config"])
                    params += ["--variant", vrn_file]
                # For large inputs, reduce memory usage by batching
                # https://github.com/bcbio/bcbio-nextgen/issues/2852
                if len(vrn_files) > 200:
                    params += ["--batch-size", "50"]
                memscale = {"magnitude": 0.9 * cores, "direction": "increase"} if cores > 1 else None
                broad_runner.run_gatk(params, memscale=memscale)
    return out_dir

def _incomplete_genomicsdb(dbdir):
    """Check if a GenomicsDB output is incomplete and we should regenerate.

    Works around current inability to move GenomicsDB outputs and support
    transactional directories.
    """
    for test_file in ["callset.json", "vidmap.json", "genomicsdb_array/genomicsdb_meta.json"]:
        if not os.path.exists(os.path.join(dbdir, test_file)):
            return True
    return False

def _run_genotype_gvcfs_genomicsdb(genomics_db, region, out_file, data):
    """GenotypeGVCFs from a merged GenomicsDB input: GATK4.
            ropts += [str(x) for x in resources.get("options", [])]

    No core scaling -- not yet supported in GATK4.
    """
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            broad_runner = broad.runner_from_config(data["config"])
            params = ["-T", "GenotypeGVCFs",
                      "--variant", "gendb://%s" % genomics_db,
                      "-R", dd.get_ref_file(data),
                      "--output", tx_out_file,
                      "-L", bamprep.region_to_gatk(region)]
            params += ["-ploidy", str(ploidy.get_ploidy([data], region))]
            # Avoid slow genotyping runtimes with improved quality score calculation in GATK4
            # https://gatkforums.broadinstitute.org/gatk/discussion/11471/performance-troubleshooting-tips-for-genotypegvcfs/p1
            params += ["--use-new-qual-calculator"]
            resources = config_utils.get_resources("gatk", data["config"])
            params += [str(x) for x in resources.get("options", [])]
            cores = dd.get_cores(data)
            memscale = {"magnitude": 0.9 * cores, "direction": "increase"} if cores > 1 else None
            broad_runner.run_gatk(params, memscale=memscale)
    return vcfutils.bgzip_and_index(out_file, data["config"])

# ## gVCF joint genotype calling -- GATK3

def _run_genotype_gvcfs_gatk3(data, region, vrn_files, ref_file, out_file):
    """Performs genotyping of gVCFs into final VCF files.
    """
    if not utils.file_exists(out_file):
        broad_runner = broad.runner_from_config(data["config"])
        with file_transaction(data, out_file) as tx_out_file:
            assoc_files = tz.get_in(("genome_resources", "variation"), data, {})
            if not assoc_files: assoc_files = {}
            params = ["-T", "GenotypeGVCFs",
                      "-R", ref_file, "-o", tx_out_file,
                      "-L", bamprep.region_to_gatk(region),
                      "--max_alternate_alleles", "4"]
            for vrn_file in vrn_files:
                params += ["--variant", vrn_file]
            if assoc_files.get("dbsnp"):
                params += ["--dbsnp", assoc_files["dbsnp"]]
            broad_runner.new_resources("gatk-haplotype")
            cores = dd.get_cores(data)
            if cores > 1:
                # GATK performs poorly with memory usage when parallelizing
                # with a large number of cores but makes use of extra memory,
                # so we cap at 6 cores.
                # See issue #1565 for discussion
                # Recent GATK 3.x versions also have race conditions with multiple
                # threads, so limit to 1 and keep memory available
                # https://gatkforums.broadinstitute.org/wdl/discussion/8718/concurrentmodificationexception-in-gatk-3-7-genotypegvcfs
                # params += ["-nt", str(min(6, cores))]
                memscale = {"magnitude": 0.9 * cores, "direction": "increase"}
            else:
                memscale = None
            broad_runner.run_gatk(params, memscale=memscale, parallel_gc=True)
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
            broad_runner.run_gatk(params, memscale=memscale, parallel_gc=True)
    return vcfutils.bgzip_and_index(out_file, data["config"])
