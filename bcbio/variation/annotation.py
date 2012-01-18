"""Annotated variant VCF files with additional information.

- GATK variant annotation with snpEff predicted effects.
"""
import os

from bcbio import broad
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction

# ## snpEff annotation

def annotate_effects(orig_file, snpeff_file, genome_file, config):
    """Annotate predicted variant effects using snpEff.
    """
    broad_runner = broad.runner_from_config(config)
    out_file = "%s-annotated%s" % os.path.splitext(orig_file)
    # Avoid generalization since 2.0.3 is not working
    #snpeff_file = _general_snpeff_version(snpeff_file)
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "VariantAnnotator",
                      "-R", genome_file,
                      "-A", "SnpEff",
                      "--variant", orig_file,
                      "--snpEffFile", snpeff_file,
                      "--out", tx_out_file]
            broad_runner.run_gatk(params)
    return out_file

def _fix_snpeff_version_line(line, supported_versions):
    """Change snpEff versions to supported minor releases.

    ##SnpEffVersion="2.0.3 (build 2011-10-08), by Pablo Cingolani"
    """
    start, rest = line.split('"', 1)
    version, end = rest.split(" ", 1)
    version_base = version.rsplit(".", 1)[0]
    for sv in supported_versions:
        sv_base = sv.rsplit(".", 1)[0]
        if sv_base == version_base:
            version = sv
            break
    return '%s"%s %s' % (start, version, end)

def _general_snpeff_version(snpeff_file):
    """GATK wants exact snpEff versions; allow related minor releases.
    """
    gatk_versions = ["2.0.2"]
    safe_snpeff = "%s-safev%s" % os.path.splitext(snpeff_file)
    if not file_exists(safe_snpeff):
        with file_transaction(safe_snpeff) as tx_safe:
            with open(snpeff_file) as in_handle:
                with open(safe_snpeff, "w") as out_handle:
                    for line in in_handle:
                        if line.startswith("##SnpEffVersion"):
                            line = _fix_snpeff_version_line(line, gatk_versions)
                        out_handle.write(line)
    return safe_snpeff

def annotate_nongatk_vcf(orig_file, dbsnp_file, ref_file, config):
    """Annotate a VCF file with dbSNP and standard GATK called annotations.
    """
    broad_runner = broad.runner_from_config(config)
    out_file = "%s-gatkann%s" % os.path.splitext(orig_file)
    annotations = ["BaseQualityRankSumTest", "DepthOfCoverage", "FisherStrand",
                   "GCContent", "HaplotypeScore", "HomopolymerRun",
                   "MappingQualityRankSumTest", "MappingQualityZero",
                   "QualByDepth", "ReadPosRankSumTest", "RMSMappingQuality"]
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "VariantAnnotator",
                      "-R", ref_file,
                      "--variant", orig_file,
                      "--dbsnp", dbsnp_file,
                      "--out", tx_out_file]
            for x in annotations:
                params += ["-A", x]
            broad_runner.run_gatk(params)
    return out_file
