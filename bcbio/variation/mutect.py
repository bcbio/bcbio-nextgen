"""Provide support for MuTect and other paired analysis tools."""

from distutils.version import LooseVersion
import os

from bcbio import bam, broad
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.variation.realign import has_aligned_reads
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.variation import bamprep, vcfutils
from bcbio.log import logger

_PASS_EXCEPTIONS = set(["java.lang.RuntimeException: "
                        "java.lang.IllegalArgumentException: "
                        "Comparison method violates its general contract!",
                        "java.lang.IllegalArgumentException: "
                        "Comparison method violates its general contract!"])

def _check_mutect_version(broad_runner):
    mutect_version = broad_runner.get_mutect_version()
    try:
        assert mutect_version is not None
    except AssertionError:
        logger.warn("WARNING")
        logger.warn("MuTect version could not be determined from jar file. "
                    "Please ensure you are using at least version 1.1.5, "
                    "as versions 1.1.4 and lower have known issues.")
        logger.warn("Proceeding but assuming correct version 1.1.5.")
    else:
        try:
            assert LooseVersion(mutect_version) >= LooseVersion("1.1.5")
        except AssertionError:
            message = ("MuTect 1.1.4 and lower is known to have incompatibilities "
                       "with Java < 7, and this may lead to problems in analyses. "
                       "Please use MuTect 1.1.5 or higher (note that it requires "
                       "Java 7).")
            raise ValueError(message)

def _config_params(base_config, assoc_files, region, out_file):
    """Add parameters based on configuration variables, associated files and genomic regions.
    """
    params = []
    contamination = base_config["algorithm"].get("fraction_contamination", 0)
    params += ["--fraction_contamination", contamination]
    dbsnp = assoc_files["dbsnp"]
    if dbsnp:
        params += ["--dbsnp", dbsnp]
    cosmic = assoc_files.get("cosmic")
    if cosmic:
        params += ["--cosmic", cosmic]
    variant_regions = base_config["algorithm"].get("variant_regions")
    region = subset_variant_regions(variant_regions, region, out_file)
    if region:
        params += ["-L", bamprep.region_to_gatk(region), "--interval_set_rule",
                   "INTERSECTION"]
    return params

def _mutect_call_prep(align_bams, items, ref_file, assoc_files,
                       region=None, out_file=None):
    """Preparation work for MuTect.
    """
    base_config = items[0]["config"]
    broad_runner = broad.runner_from_config(base_config, "mutect")
    _check_mutect_version(broad_runner)

    broad_runner.run_fn("picard_index_ref", ref_file)
    for x in align_bams:
        bam.index(x, base_config)

    params = ["-R", ref_file, "-T", "MuTect"]
    paired = vcfutils.get_paired_bams(align_bams, items)
    params += ["-I:tumor", paired.tumor_bam]
    params += ["--tumor_sample_name", paired.tumor_name]
    if paired.normal_bam is not None:
        params += ["-I:normal", paired.normal_bam]
        params += ["--normal_sample_name", paired.normal_name]
    if paired.normal_panel is not None:
        params += ["--normal_panel", paired.normal_panel]
    params += _config_params(base_config, assoc_files, region, out_file)
    return broad_runner, params

def mutect_caller(align_bams, items, ref_file, assoc_files, region=None,
                  out_file=None):
    """Run the MuTect paired analysis algorithm.
    """
    if out_file is None:
        out_file = "%s-paired-variants.vcf" % os.path.splitext(align_bams[0])[0]
    if not file_exists(out_file):
        broad_runner, params = \
            _mutect_call_prep(align_bams, items, ref_file, assoc_files,
                                   region, out_file)
        if (not isinstance(region, (list, tuple)) and
              not all(has_aligned_reads(x, region) for x in align_bams)):
                vcfutils.write_empty_vcf(out_file)
                return
        with file_transaction(out_file) as tx_out_file:
            # Rationale: MuTect writes another table to stdout, which we don't need
            params += ["--vcf", tx_out_file, "-o", os.devnull]
            broad_runner.run_mutect(params)
    return out_file
