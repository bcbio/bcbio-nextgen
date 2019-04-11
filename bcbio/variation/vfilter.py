"""Cutoff-based soft filtering of genomic variants.
"""
from distutils.version import LooseVersion
import math
import os
import shutil

import numpy
import toolz as tz
import yaml

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do, programs
from bcbio.variation import vcfutils

# ## General functionality

def cutoff_w_expression(vcf_file, expression, data, name="+", filterext="",
                      extra_cmd="", limit_regions="variant_regions"):
    """Perform cutoff-based soft filtering using bcftools expressions like %QUAL < 20 || DP < 4.
    """
    base, ext = utils.splitext_plus(vcf_file)
    out_file = "{base}-filter{filterext}{ext}".format(**locals())
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            if vcfutils.vcf_has_variants(vcf_file):
                bcftools = config_utils.get_program("bcftools", data["config"])
                bgzip_cmd = "| bgzip -c" if out_file.endswith(".gz") else ""
                intervals = ""
                if limit_regions == "variant_regions":
                    variant_regions = dd.get_variant_regions(data)
                    if variant_regions:
                        intervals = "-T %s" % vcfutils.bgzip_and_index(variant_regions, data["config"])
                cmd = ("{bcftools} filter -O v {intervals} --soft-filter '{name}' "
                       "-e '{expression}' -m '+' {vcf_file} {extra_cmd} {bgzip_cmd} > {tx_out_file}")
                do.run(cmd.format(**locals()),
                       "Cutoff-based soft filtering %s with %s" % (vcf_file, expression), data)
            else:
                shutil.copy(vcf_file, out_file)
    if out_file.endswith(".vcf.gz"):
        out_file = vcfutils.bgzip_and_index(out_file, data["config"])
    return out_file

# ## Caller specific

def freebayes(in_file, ref_file, vrn_files, data):
    """FreeBayes filters: cutoff-based soft filtering.
    """
    out_file = _freebayes_cutoff(in_file, data)
    #out_file = _freebayes_custom(in_file, ref_file, data)
    return out_file

def _freebayes_custom(in_file, ref_file, data):
    """Custom FreeBayes filtering using bcbio.variation, tuned to human NA12878 results.

    Experimental: for testing new methods.
    """
    if vcfutils.get_paired_phenotype(data):
        return None
    config = data["config"]
    bv_ver = programs.get_version("bcbio_variation", config=config)
    if LooseVersion(bv_ver) < LooseVersion("0.1.1"):
        return None
    out_file = "%s-filter%s" % os.path.splitext(in_file)
    if not utils.file_exists(out_file):
        tmp_dir = utils.safe_makedir(os.path.join(os.path.dirname(in_file), "tmp"))
        resources = config_utils.get_resources("bcbio_variation", config)
        jvm_opts = resources.get("jvm_opts", ["-Xms750m", "-Xmx2g"])
        java_args = ["-Djava.io.tmpdir=%s" % tmp_dir]
        cmd = ["bcbio-variation"] + jvm_opts + java_args + \
              ["variant-filter", "freebayes", in_file, ref_file]
        do.run(cmd, "Custom FreeBayes filtering using bcbio.variation")
    return out_file

def _freebayes_cutoff(in_file, data):
    """Perform filtering of FreeBayes results, flagging low confidence calls.

    Filters using cutoffs on low depth based on Meynert et al's work modeling sensitivity
    of homozygote and heterozygote calling on depth:

    http://www.ncbi.nlm.nih.gov/pubmed/23773188

    and high depth heterozygote SNP filtering based on Heng Li's work
    evaluating variant calling artifacts:

    http://arxiv.org/abs/1404.0929

    Tuned based on NA12878 call comparisons to Genome in a Bottle reference genome.
    """
    if not vcfutils.vcf_has_variants(in_file):
        base, ext = utils.splitext_plus(in_file)
        out_file = "{base}-filter{ext}".format(**locals())
        if not utils.file_exists(out_file):
            shutil.copy(in_file, out_file)
        if out_file.endswith(".vcf.gz"):
            out_file = vcfutils.bgzip_and_index(out_file, data["config"])
        return out_file

    depth_thresh, qual_thresh = None, None
    if _do_high_depth_filter(data):
        stats = _calc_vcf_stats(in_file)
        if stats["avg_depth"] > 0:
            depth_thresh = int(math.ceil(stats["avg_depth"] + 3 * math.pow(stats["avg_depth"], 0.5)))
            qual_thresh = depth_thresh * 2.0  # Multiplier from default GATK QD cutoff filter
    filters = ('(AF[0] <= 0.5 && (max(FORMAT/DP) < 4 || (max(FORMAT/DP) < 13 && %QUAL < 10))) || '
               '(AF[0] > 0.5 && (max(FORMAT/DP) < 4 && %QUAL < 50))')
    if depth_thresh:
        filters += ' || (%QUAL < {qual_thresh} && max(FORMAT/DP) > {depth_thresh} && AF[0] <= 0.5)'.format(**locals())
    return cutoff_w_expression(in_file, filters, data, name="FBQualDepth")

def _do_high_depth_filter(data):
    """Check if we should do high depth filtering -- only on germline non-regional calls.
    """
    is_genome = tz.get_in(["config", "algorithm", "coverage_interval"], data, "").lower() == "genome"
    is_paired = vcfutils.get_paired_phenotype(data)
    return is_genome and not is_paired

def _calc_vcf_stats(in_file):
    """Calculate statistics on VCF for filtering, saving to a file for quick re-runs.
    """
    out_file = "%s-stats.yaml" % utils.splitext_plus(in_file)[0]
    if not utils.file_exists(out_file):
        stats = {"avg_depth": _average_called_depth(in_file)}
        with open(out_file, "w") as out_handle:
            yaml.safe_dump(stats, out_handle, default_flow_style=False, allow_unicode=False)
        return stats
    else:
        with open(out_file) as in_handle:
            stats = yaml.safe_load(in_handle)
        return stats

def _average_called_depth(in_file):
    """Retrieve the average depth of called reads in the provided VCF.
    """
    import cyvcf2
    depths = []
    for rec in cyvcf2.VCF(str(in_file)):
        d = rec.INFO.get("DP")
        if d is not None:
            depths.append(int(d))
    if len(depths) > 0:
        return int(math.ceil(numpy.mean(depths)))
    else:
        return 0

def platypus(in_file, data):
    """Filter Platypus calls, removing Q20 filter and replacing with depth and quality based filter.

    Platypus uses its own VCF nomenclature: TC == DP, FR == AF

    Platypus gVCF output appears to have an 0/1 index problem so the reference block
    regions are 1 base outside regions of interest. We avoid limiting regions during
    filtering when using it.
    """
    filters = ('(FR[0] <= 0.5 && TC < 4 && %QUAL < 20) || '
               '(TC < 13 && %QUAL < 10) || '
               '(FR[0] > 0.5 && TC < 4 && %QUAL < 50)')
    limit_regions = "variant_regions" if not vcfutils.is_gvcf_file(in_file) else None
    return cutoff_w_expression(in_file, filters, data, name="PlatQualDepth",
                               extra_cmd="| sed 's/\\tQ20\\t/\\tPASS\\t/'", limit_regions=limit_regions)

def samtools(in_file, data):
    """Filter samtools calls based on depth and quality, using similar approaches to FreeBayes.
    """
    filters = ('((AC[0] / AN) <= 0.5 && max(FORMAT/DP) < 4 && %QUAL < 20) || '
               '(max(FORMAT/DP) < 13 && %QUAL < 10) || '
               '((AC[0] / AN) > 0.5 && max(format/DP) < 4 && %QUAL < 50)')
    return cutoff_w_expression(in_file, filters, data, name="stQualDepth")

def _gatk_general():
    """General filters useful for both GATK SNPs and indels.

    Remove low quality, low allele fraction variants at the ends of reads.
    Generally useful metric identified by looking at 10x data.
    https://community.10xgenomics.com/t5/Genome-Exome-Forum/Best-practices-for-trimming-adapters-when-variant-calling/m-p/473
    https://github.com/bcbio/bcbio_validations/tree/master/gatk4#10x-adapter-trimming--low-frequency-allele-filter
    """
    return ["(QD < 10.0 && AD[0:1] / (AD[0:1] + AD[0:0]) < 0.25 && ReadPosRankSum < 0.0)"]

def gatk_snp_cutoff(in_file, data):
    """Perform cutoff-based soft filtering on GATK SNPs using best-practice recommendations.

    We have a more lenient mapping quality (MQ) filter compared to GATK defaults.
    The recommended filter (MQ < 40) is too stringent, so we adjust to 30:
    http://imgur.com/a/oHRVB

    QD and FS are not calculated when generating gVCF output:
    https://github.com/broadgsa/gatk-protected/blob/e91472ddc7d58ace52db0cab4d70a072a918d64c/protected/gatk-tools-protected/src/main/java/org/broadinstitute/gatk/tools/walkers/haplotypecaller/HaplotypeCaller.java#L300

    The extra command removes escaped quotes in the VCF output which
    pyVCF fails on.

    Does not use the GATK best practice recommend SOR filter (SOR > 3.0) as it
    has a negative impact on sensitivity relative to precision:

    https://github.com/bcbio/bcbio_validations/tree/master/gatk4#na12878-hg38
    """
    filters = ["MQRankSum < -12.5", "ReadPosRankSum < -8.0"]
    # GATK Haplotype caller (v2.2) appears to have much larger HaplotypeScores
    # resulting in excessive filtering, so avoid this metric
    variantcaller = utils.get_in(data, ("config", "algorithm", "variantcaller"))
    if variantcaller not in ["gatk-haplotype", "haplotyper"]:
        filters.append("HaplotypeScore > 13.0")
    # Additional filter metrics, unless using raw GATK HaplotypeCaller or Sentieon gVCFs
    if not (vcfutils.is_gvcf_file(in_file) and variantcaller in ["gatk-haplotype", "haplotyper"]):
        filters += ["QD < 2.0"]
        filters += ["FS > 60.0"]
        filters += _gatk_general()
        filters += ["MQ < 30.0"]
    return cutoff_w_expression(in_file, 'TYPE="snp" && (%s)' % " || ".join(filters), data, "GATKCutoffSNP", "SNP",
                               extra_cmd=r"""| sed 's/\\"//g'""")

def gatk_indel_cutoff(in_file, data):
    """Perform cutoff-based soft filtering on GATK indels using best-practice recommendations.
    """
    filters = ["ReadPosRankSum < -20.0"]
    variantcaller = utils.get_in(data, ("config", "algorithm", "variantcaller"))
    # Additional filter metrics, unless using raw GATK HaplotypeCaller or Sentieon gVCFs
    if not (vcfutils.is_gvcf_file(in_file) and variantcaller in ["gatk-haplotype", "haplotyper"]):
        filters += ["QD < 2.0"]
        filters += ["FS > 200.0"]
        filters += ["SOR > 10.0"]
        filters += _gatk_general()
    return cutoff_w_expression(in_file, 'TYPE="indel" && (%s)' % " || ".join(filters), data, "GATKCutoffIndel",
                               "INDEL", extra_cmd=r"""| sed 's/\\"//g'""")
