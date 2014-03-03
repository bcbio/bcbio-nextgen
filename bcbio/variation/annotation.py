"""Annotated variant VCF files with additional information.

- GATK variant annotation with snpEff predicted effects.
"""
import os

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.variation import vcfutils

def get_gatk_annotations(config):
    broad_runner = broad.runner_from_config(config)
    anns = ["BaseQualityRankSumTest", "FisherStrand",
            "GCContent", "HaplotypeScore", "HomopolymerRun",
            "MappingQualityRankSumTest", "MappingQualityZero",
            "QualByDepth", "ReadPosRankSumTest", "RMSMappingQuality",
            "DepthPerAlleleBySample"]
    if broad_runner.gatk_type() == "restricted":
        anns += ["Coverage"]
    else:
        anns += ["DepthOfCoverage"]
    return anns

def annotate_nongatk_vcf(orig_file, bam_files, dbsnp_file, ref_file, config):
    """Annotate a VCF file with dbSNP and standard GATK called annotations.
    """
    orig_file = vcfutils.bgzip_and_index(orig_file, config)
    out_file = "%s-gatkann%s" % utils.splitext_plus(orig_file)
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            # Avoid issues with incorrectly created empty GATK index files.
            # Occurs when GATK cannot lock shared dbSNP database on previous run
            idx_file = orig_file + ".idx"
            if os.path.exists(idx_file) and not utils.file_exists(idx_file):
                os.remove(idx_file)
            annotations = get_gatk_annotations(config)
            params = ["-T", "VariantAnnotator",
                      "-R", ref_file,
                      "--variant", orig_file,
                      "--dbsnp", dbsnp_file,
                      "--out", tx_out_file,
                      "-L", orig_file]
            for bam_file in bam_files:
                params += ["-I", bam_file]
            for x in annotations:
                params += ["-A", x]
            broad_runner = broad.runner_from_config(config)
            broad_runner.run_gatk(params, memory_retry=True)
    vcfutils.bgzip_and_index(out_file, config)
    return out_file
