"""GATK variant calling -- MuTect2.
"""
from distutils.version import LooseVersion

import toolz as tz

from bcbio import bam, broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.pipeline import datadict as dd
from bcbio.variation import annotation, bamprep, vcfutils, ploidy
from bcbio.variation.vcfutils import bgzip_and_index

def _shared_gatk_call_prep(align_bams, items, ref_file, dbsnp, cosmic, region, out_file):
    """Shared preparation work for GATK variant calling.
    """
    data = items[0]
    config = data["config"]
    broad_runner = broad.runner_from_path("picard", config)
    broad_runner.run_fn("picard_index_ref", ref_file)
    for x in align_bams:
        bam.index(x, config)
    params = ["-R", ref_file]
    coverage_depth_min = tz.get_in(["algorithm", "coverage_depth_min"], config)
    if coverage_depth_min and coverage_depth_min < 4:
        confidence = "4.0"
        params += ["--standard_min_confidence_threshold_for_calling", confidence,
                   "--standard_min_confidence_threshold_for_emitting", confidence]
    for a in annotation.get_gatk_annotations(config):
        params += ["--annotation", a]
    for x in align_bams:
        bam.index(x, config)

    paired = vcfutils.get_paired_bams(align_bams, items)
    if not paired:
        raise ValueError("Specified MuTect calling but 'tumor' phenotype not present in batch\n"
                         "https://bcbio-nextgen.readthedocs.org/en/latest/contents/"
                         "pipelines.html#cancer-variant-calling\n"
                         "for samples: %s" % ", " .join([dd.get_sample_name(x) for x in items]))
    params += ["-I:tumor", paired.tumor_bam]
    if paired.normal_bam is not None:
        params += ["-I:normal", paired.normal_bam]
    if paired.normal_panel is not None:
        params += ["--normal_panel", paired.normal_panel]
    if dbsnp:
        params += ["--dbsnp", dbsnp]
    if cosmic:
        params += ["--cosmic", cosmic]
    variant_regions = tz.get_in(["algorithm", "variant_regions"], config)
    region = subset_variant_regions(variant_regions, region, out_file, items)
    if region:
        params += ["-L", bamprep.region_to_gatk(region), "--interval_set_rule", "INTERSECTION"]
    broad_runner = broad.runner_from_config(config)
    return broad_runner, params
   
def mutect2_caller(align_bams, items, ref_file, assoc_files,
                       region=None, out_file=None):
    """Call variation with GATK's MuTect2.

    This requires the full non open-source version of GATK 3.5+.
    """
    if out_file is None:
        out_file = "%s-variants.vcf.gz" % utils.splitext_plus(align_bams[0])[0]
    if not utils.file_exists(out_file):
        broad_runner, params = \
            _shared_gatk_call_prep(align_bams, items,
                                   ref_file, assoc_files.get("dbsnp"), assoc_files.get("cosmic"),
                                   region, out_file)
        assert LooseVersion(broad_runner.gatk_major_version()) >= LooseVersion("3.5"), \
            "Require full version of GATK 3.5+ for mutect2 calling"
        with file_transaction(items[0], out_file) as tx_out_file:
            params += ["-T", "MuTect2",
                       "-o", tx_out_file,
                       "--annotation", "ClippingRankSumTest",
                       "--annotation", "DepthPerSampleHC"]
            resources = config_utils.get_resources("mutect2", items[0]["config"])
            if "options" in resources:
                params += [str(x) for x in resources.get("options", [])]
            broad_runner.new_resources("mutect2")
            broad_runner.run_gatk(params)
    return out_file

