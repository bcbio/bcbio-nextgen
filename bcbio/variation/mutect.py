"""Provide support for MuTect and other paired analysis tools."""

from distutils.version import LooseVersion
import os
#import itertools
#from subprocess import CalledProcessError

from bcbio import bam, broad
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.variation.realign import has_aligned_reads
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.variation import bamprep, vcfutils
from bcbio.variation.vcfutils import get_paired_bams, is_paired_analysis
from bcbio.log import logger


_PASS_EXCEPTIONS = set(["java.lang.RuntimeException: "
                        "java.lang.IllegalArgumentException: "
                        "Comparison method violates its general contract!",
                        "java.lang.IllegalArgumentException: "
                        "Comparison method violates its general contract!"])


def _mutect_call_prep(align_bams, items, ref_file, assoc_files,
                       region=None, out_file=None):
    """
    Preparation work for MuTect.
    """

    #FIXME: We assume all other bits in the config are shared

    base_config = items[0]["config"]
    dbsnp = assoc_files["dbsnp"]
    cosmic = assoc_files.get("cosmic")

    broad_runner = broad.runner_from_config(base_config, "mutect")

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
            message =  ("MuTect 1.1.4 and lower is known to have incompatibilities "
                        "with Java < 7, and this may lead to problems in analyses. "
                        "Please use MuTect 1.1.5 or higher (note that it requires "
                        "Java 7).")
            raise ValueError(message)

    broad_runner.run_fn("picard_index_ref", ref_file)
    for x in align_bams:
        bam.index(x, base_config)

    variant_regions = base_config["algorithm"].get("variant_regions", None)
    contamination = base_config["algorithm"].get("fraction_contamination", 0)
    region = subset_variant_regions(variant_regions, region, out_file)

    #FIXME: Add more parameters like fraction contamination etc

    params = ["-R", ref_file, "-T", "MuTect"]
    params += ["--dbsnp", dbsnp]
    
    # mutect supports single-normal bams only (Feb2014) and multiple tumor bams
    # are concatenated prior to calling and thus there are only two samples 
    # ever in the vcf output
    # 
    # check if a tumor/normal pair is defined or whether it's tumor only 
    if is_paired_analysis(align_bams, items):
        paired = get_paired_bams(align_bams, items)
        params += ["-I:normal", paired.normal_bam]
        params += ["--normal_sample_name", paired.normal_sample_name]
        tumor_bam = paired.tumor_bam
        tumor_sample_name = paired.tumor_sample_name
        panel_of_normals_vcf = paired.normal_panel
    else:
        # sanity check - there should be no singleton normals or more than 2 bams
        if (len(align_bams)>1 or items[0]["metadata"]["phenotype"] == "normal"):
            message =  ("MuTect currently supports tumor/normal and tumor only modes."
                        "Provide one sample with tumor phenotype or two samples with "
                        "one tumor and one normal phenotype. ")
            raise ValueError(message)
        tumor_bam = align_bams[0]
        tumor_sample_name = items[0]["name"][1]
        panel_of_normals_vcf = items[0]["metadata"].get("normal_panel",None)

    params += ["-I:tumor", tumor_bam]
    params += ["--tumor_sample_name", tumor_sample_name]
    params += ["--fraction_contamination", contamination]
    
    if panel_of_normals_vcf is not None:
        params += ["--normal_panel", panel_of_normals_vcf]
    
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

            broad_runner.run_mutect(params)

    return out_file
