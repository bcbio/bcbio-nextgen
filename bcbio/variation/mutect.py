"""Provide support for MuTect and other paired analysis tools."""

import os
import itertools
from subprocess import CalledProcessError

from bcbio import broad
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.variation.realign import has_aligned_reads
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.variation import bamprep, vcfutils

_PASS_EXCEPTIONS = set(["java.lang.RuntimeException: "
                        "java.lang.IllegalArgumentException: "
                        "Comparison method violates its general contract!",
                        "java.lang.IllegalArgumentException: "
                        "Comparison method violates its general contract!"])


def _parse_gatk_java_error_string(error_string):

    """Parse the GATK error string to get the stack trace"""

    for line in error_string.split("##### ERROR "):
        line = line.strip()
        if "stack trace" in line:
            line = line.split("\n")
            # The name of the exception is immediately after "stack trace"
            java_error = line[1].strip()

            return java_error


def _mutect_call_prep(align_bams, items, ref_file, assoc_files,
                       region=None, out_file=None):
    """
    Preparation work for MuTect.
    """

    #FIXME: We assume all other bits in the config are shared

    base_config = items[0]["config"]
    dbsnp = assoc_files.dbsnp
    cosmic = None

    if hasattr(assoc_files, "cosmic"):
        cosmic = assoc_files.cosmic

    broad_runner = broad.runner_from_config(base_config, "mutect")

    broad_runner.run_fn("picard_index_ref", ref_file)
    for x in align_bams:
        broad_runner.run_fn("picard_index", x)

    variant_regions = base_config["algorithm"].get("variant_regions", None)
    contamination = base_config["algorithm"].get("fraction_contamination", 0)
    region = subset_variant_regions(variant_regions, region, out_file)

    #FIXME: Add more parameters like fraction contamination etc

    params = ["-R", ref_file, "-T", "MuTect"]
    params += ["--dbsnp", dbsnp]

    tumor_bam = None
    normal_bam = None

    for bamfile, item in itertools.izip(align_bams, items):

        metadata = item["metadata"]

        if metadata["phenotype"] == "normal":
            normal_bam = bamfile
            normal_sample_name = item["name"][1]
        elif metadata["phenotype"] == "tumor":
            tumor_bam = bamfile
            tumor_sample_name = item["name"][1]

    if tumor_bam is None or normal_bam is None:
        raise ValueError("Missing phenotype definition (tumor or normal) "
                         "in samples")

    params += ["-I:normal", normal_bam]
    params += ["-I:tumor", tumor_bam]
    params += ["--tumor_sample_name", tumor_sample_name]
    params += ["--normal_sample_name", normal_sample_name]
    params += ["--fraction_contamination", contamination]

    if cosmic is not None:
        params += ["--cosmic", cosmic]

    if region:
        params += ["-L", bamprep.region_to_gatk(region), "--interval_set_rule",
                   "INTERSECTION"]

    return broad_runner, params


def mutect_caller(align_bams, items, ref_file, assoc_files, region=None,
                  out_file=None):

    """Run the MuTect paired analysis algorithm."""

    if out_file is None:
        out_file = "%s-paired-variants.vcf" % os.path.splitext(
            align_bams[0])[0]

    if not file_exists(out_file):
        broad_runner, params = \
            _mutect_call_prep(align_bams, items, ref_file, assoc_files,
                                   region, out_file)

        if (not isinstance(region, (list, tuple)) and
            not all(has_aligned_reads(x, region) for x in align_bams)):

                vcfutils.write_empty_vcf(out_file)
                return

        with file_transaction(out_file) as tx_out_file:
            # Rationale: MuTect writes another table to stdout,
            # which we don't need
            params += ["--vcf", tx_out_file, "-o", os.devnull]
            try:
                broad_runner.run_mutect(params)
            except CalledProcessError as error:
                java_exception = _parse_gatk_java_error_string(error.cmd)
                #HACK: Currently MuTect bails out on certain small BAM files
                # Until the issue is fixed by Broad, this specific exception
                # will be ignored. All the other exceptions will be raised
                # correctly.
                if java_exception in _PASS_EXCEPTIONS:
                    vcfutils.write_empty_vcf(tx_out_file)
                    return
                else:
                    raise

    return out_file
