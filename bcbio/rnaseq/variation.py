import os
from bcbio import utils
from bcbio.utils import file_exists, get_R_exports, safe_makedir
import bcbio.pipeline.datadict as dd
from bcbio.pipeline import config_utils
from bcbio.ngsalign.postalign import dedup_bam
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.variation import vardict
from bcbio import broad, bam
from bcbio.variation import gatk, vcfutils
from bcbio.rnaseq import gtf

def rnaseq_gatk_variant_calling(data):
    data = dd.set_deduped_bam(data, dedup_bam(dd.get_work_bam(data), data))
    data = gatk_splitreads(data)
    data = gatk_rnaseq_calling(data)
    return data

def gatk_splitreads(data):
    """
    use GATK to split reads with Ns in the CIGAR string, hard clipping regions
    that end up in introns
    """
    broad_runner = broad.runner_from_config(dd.get_config(data))
    ref_file = dd.get_ref_file(data)
    deduped_bam = dd.get_deduped_bam(data)
    base, ext = os.path.splitext(deduped_bam)
    split_bam = base + ".splitN" + ext
    if file_exists(split_bam):
        data = dd.set_split_bam(data, split_bam)
        return data
    gatk_type = broad_runner.gatk_type()
    with file_transaction(data, split_bam) as tx_split_bam:
        params = ["-T", "SplitNCigarReads",
                  "-R", ref_file,
                  "-I", deduped_bam]
        if gatk_type == "gatk4":
            params += ["--output", tx_split_bam]
        else:
            params += ["-rf", "ReassignOneMappingQuality",
                       "-RMQF", "255",
                       "-RMQT", "60",
                       "-rf", "UnmappedRead",
                       "-U", "ALLOW_N_CIGAR_READS",
                       "-o", tx_split_bam]
            if dd.get_quality_format(data) == "illumina":
                params += ["--fix_misencoded_quality_scores", "-fixMisencodedQuals"]
        broad_runner.run_gatk(params)
    bam.index(split_bam, dd.get_config(data))
    data = dd.set_split_bam(data, split_bam)
    return data

def gatk_rnaseq_calling(data):
    """Use GATK to perform gVCF variant calling on RNA-seq data
    """
    data = utils.deepish_copy(data)
    tools_on = dd.get_tools_on(data)
    if not tools_on:
        tools_on = []
    tools_on.append("gvcf")
    data = dd.set_tools_on(data, tools_on)
    data = dd.set_jointcaller(data, ["%s-joint" % v for v in dd.get_variantcaller(data)])
    out_file = os.path.join(utils.safe_makedir(os.path.join("variation", "rnaseq", "gatk-haplotype")),
                            "%s-gatk-haplotype.vcf.gz" % dd.get_sample_name(data))
    out_file = gatk.haplotype_caller([dd.get_split_bam(data)], [data], dd.get_ref_file(data), {},
                                     out_file=out_file)
    return dd.set_vrn_file(data, out_file)

def rnaseq_vardict_variant_calling(data):
    sample = dd.get_sample_name(data)
    variation_dir = os.path.join(dd.get_work_dir(data), "variation")
    safe_makedir(variation_dir)
    out_file = os.path.join(variation_dir, sample + "-vardict.vcf.gz")
    if file_exists(out_file):
        data = dd.set_vrn_file(data, out_file)
        return data
    vardict_cmd = vardict.get_vardict_command(data)
    strandbias = "teststrandbias.R"
    var2vcf = "var2vcf_valid.pl"
    vcfstreamsort = config_utils.get_program("vcfstreamsort", data)
    compress_cmd = "| bgzip -c"
    freq = float(dd.get_min_allele_fraction(data, 20) / 100.0)
    var2vcf_opts = "-v 50"
    fix_ambig = vcfutils.fix_ambiguous_cl()
    remove_dup = vcfutils.remove_dup_cl()
    r_setup = get_R_exports()
    ref_file = dd.get_ref_file(data)
    bamfile = dd.get_work_bam(data)
    bed_file = gtf.gtf_to_bed(dd.get_gtf_file(data))
    opts = " -c 1 -S 2 -E 3 -g 4 "
    resources = config_utils.get_resources("vardict", data)
    if resources.get("options"):
        opts += " ".join([str(x) for x in resources["options"]])
    with file_transaction(data, out_file) as tx_out_file:
        jvm_opts = vardict._get_jvm_opts(data, tx_out_file)
        cmd = ("{r_setup} && {jvm_opts}{vardict_cmd} -G {ref_file} -f {freq} "
                "-N {sample} -b {bamfile} {opts} {bed_file} "
                "| {strandbias}"
                "| {var2vcf} -N {sample} -E -f {freq} {var2vcf_opts} "
                "| {fix_ambig} | {remove_dup} | {vcfstreamsort} {compress_cmd} "
                "> {tx_out_file}")
        message = "Calling RNA-seq variants with VarDict"
        do.run(cmd.format(**locals()), message)
    data = dd.set_vrn_file(data, out_file)
    return data

def gatk_filter_rnaseq(vrn_file, data):
    """
    this incorporates filters listed here, dropping clusters of variants
    within a 35 nucleotide window, high fischer strand values and low
    quality by depth
    https://software.broadinstitute.org/gatk/guide/article?id=3891
    java -jar GenomeAnalysisTK.jar -T VariantFiltration -R hg_19.fasta -V
    input.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0"
    -filterName QD -filter "QD < 2.0" -o output.vcf
    """
    out_file = "%s-filter%s" % utils.splitext_plus(vrn_file)
    if not file_exists(out_file):
        ref_file = dd.get_ref_file(data)
        with file_transaction(data, out_file) as tx_out_file:
            params = ["-T", "VariantFiltration",
                      "-R", ref_file,
                      "-V", vrn_file,
                      "--clusterWindowSize", "35",
                      "--clusterSize", "3",
                      "--filterExpression", "\"'FS > 30.0'\"",
                      "--filterName", "FS",
                      "--filterExpression", "\"'QD < 2.0'\"",
                      "--filterName", "QD",
                      "-o", tx_out_file]
            jvm_opts = broad.get_gatk_framework_opts(dd.get_config(data), os.path.dirname(tx_out_file))
            do.run(broad.gatk_cmd("gatk-framework", jvm_opts, params), "Filter RNA-seq variants.")
    return out_file
