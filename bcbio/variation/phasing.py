"""Approaches for calculating haplotype phasing of variants.
"""
import os

from bcbio import broad
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import shared
from bcbio.variation import bamprep

def has_variants(vcf_file):
    with open(vcf_file) as in_handle:
        for line in in_handle:
            if not line.startswith("#"):
                return True
    return False

def read_backed_phasing(vcf_file, bam_files, genome_file, region, config):
    """Phase variants using GATK's read-backed phasing.
    http://www.broadinstitute.org/gatk/gatkdocs/
    org_broadinstitute_sting_gatk_walkers_phasing_ReadBackedPhasing.html
    """
    if has_variants(vcf_file):
        broad_runner = broad.runner_from_config(config)
        out_file = "%s-phased%s" % os.path.splitext(vcf_file)
        if not file_exists(out_file):
            with file_transaction(config, out_file) as tx_out_file:
                params = ["-T", "ReadBackedPhasing",
                          "-R", genome_file,
                          "--variant", vcf_file,
                          "--out", tx_out_file,
                          "--downsample_to_coverage", "250",
                          "--downsampling_type", "BY_SAMPLE"]
                for bam_file in bam_files:
                    params += ["-I", bam_file]
                variant_regions = config["algorithm"].get("variant_regions", None)
                region = shared.subset_variant_regions(variant_regions, region, out_file)
                if region:
                    params += ["-L", bamprep.region_to_gatk(region),
                               "--interval_set_rule", "INTERSECTION"]
                broad_runner.run_gatk(params)
        return out_file
    else:
        return vcf_file
