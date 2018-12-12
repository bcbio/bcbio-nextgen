import os
from bcbio import utils
from bcbio.utils import file_exists, get_R_exports
from bcbio.bam import ref
from bcbio.heterogeneity import chromhacks
import bcbio.pipeline.datadict as dd
from bcbio.pipeline import config_utils, shared
from bcbio.ngsalign.postalign import dedup_bam
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.variation import vardict
from bcbio import broad, bam
from bcbio.variation import gatk, vcfutils
from bcbio.structural import regions
from bcbio.variation.bedutils import get_padded_bed_file
from bcbio.log import logger
pybedtools = utils.LazyImport("pybedtools")

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

def _setup_variant_regions(data, out_dir):
    """Ensure we have variant regions for calling, using transcript if not present.

    Respects noalt_calling by removing additional contigs to improve
    speeds.
    """
    vr_file = dd.get_variant_regions(data)
    if not vr_file:
        vr_file = regions.get_sv_bed(data, "transcripts", out_dir=out_dir)
    contigs = set([c.name for c in ref.file_contigs(dd.get_ref_file(data))])
    out_file = os.path.join(utils.safe_makedir(os.path.join(dd.get_work_dir(data), "bedprep")),
                            "%s-rnaseq_clean.bed" % utils.splitext_plus(os.path.basename(vr_file))[0])
    if not utils.file_uptodate(out_file, vr_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                with shared.bedtools_tmpdir(data):
                    for r in pybedtools.BedTool(vr_file):
                        if r.chrom in contigs:
                            if chromhacks.is_nonalt(r.chrom):
                                out_handle.write(str(r))
    data = dd.set_variant_regions(data, out_file)
    return data

def gatk_rnaseq_calling(data):
    """Use GATK to perform gVCF variant calling on RNA-seq data
    """
    from bcbio.bam import callable
    data = utils.deepish_copy(data)
    tools_on = dd.get_tools_on(data)
    if not tools_on:
        tools_on = []
    tools_on.append("gvcf")
    data = dd.set_tools_on(data, tools_on)
    data = dd.set_jointcaller(data, ["%s-joint" % v for v in dd.get_variantcaller(data)])
    out_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data),
                                              "variation", "rnaseq", "gatk-haplotype"))
    data = _setup_variant_regions(data, out_dir)
    out_file = os.path.join(out_dir, "%s-gatk-haplotype.vcf.gz" % dd.get_sample_name(data))
    if not utils.file_exists(out_file):
        region_files = []
        regions = []
        for cur_region in callable.get_split_regions(dd.get_variant_regions(data), data):
            str_region = "_".join([str(x) for x in cur_region])
            region_file = os.path.join(utils.safe_makedir(os.path.join(dd.get_work_dir(data),
                                                                    "variation", "rnaseq", "gatk-haplotype",
                                                                    "regions")),
                                    "%s-%s-gatk-haplotype.vcf.gz" % (dd.get_sample_name(data), str_region))
            region_file = gatk.haplotype_caller([dd.get_split_bam(data)], [data], dd.get_ref_file(data), {},
                                                region=cur_region, out_file=region_file)
            region_files.append(region_file)
            regions.append(cur_region)
        out_file = vcfutils.concat_variant_files(region_files, out_file, regions,
                                                 dd.get_ref_file(data), data["config"])
    return dd.set_vrn_file(data, out_file)

def rnaseq_vardict_variant_calling(data):
    sample = dd.get_sample_name(data)
    out_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data),
                                              "variation", "rnaseq", "vardict"))
    out_file = os.path.join(out_dir, sample + "-vardict.vcf.gz")
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
    data = _setup_variant_regions(data, out_dir)
    opts, _ = vardict._vardict_options_from_config([data], data["config"], out_file, dd.get_variant_regions(data),
                                                   is_rnaseq=True)
    cores = dd.get_num_cores(data)
    if cores and cores > 1:
        opts += " -th %s" % str(cores)
    with file_transaction(data, out_file) as tx_out_file:
        jvm_opts = vardict._get_jvm_opts(data, tx_out_file)
        cmd = ("{r_setup} && {jvm_opts}{vardict_cmd} -G {ref_file} -f {freq} "
               "-N {sample} -b {bamfile} {opts} "
               "| {strandbias}"
               "| {var2vcf} -N {sample} -E -f {freq} {var2vcf_opts} "
               "| {fix_ambig} | {remove_dup} | {vcfstreamsort} {compress_cmd} "
               "> {tx_out_file}")
        message = "Calling RNA-seq variants with VarDict"
        do.run(cmd.format(**locals()), message)
    out_file = vcfutils.bgzip_and_index(out_file, data["config"])
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
            params = ["VariantFiltration",
                      "-R", ref_file,
                      "-V", vrn_file,
                      "--cluster-window-size", "35",
                      "--cluster-size", "3",
                      "--filter-expression", "'FS > 30.0'",
                      "--filter-name", "FS",
                      "--filter-expression", "'QD < 2.0'",
                      "--filter-name", "QD",
                      "--output", tx_out_file]
            # Use GATK4 for filtering, tools_off is for variant calling
            config = utils.deepish_copy(dd.get_config(data))
            if "gatk4" in dd.get_tools_off({"config": config}):
                config["algorithm"]["tools_off"].remove("gatk4")
            jvm_opts = broad.get_gatk_opts(config, os.path.dirname(tx_out_file))
            do.run(broad.gatk_cmd("gatk", jvm_opts, params, config), "Filter RNA-seq variants.")
    return out_file

def filter_junction_variants(vrn_file, data):
    """
    filter out variants within 10 basepairs of a splice junction, these are
    very prone to being false positives with RNA-seq data
    """
    SJ_BP_MASK = 10
    vrn_dir = os.path.dirname(vrn_file)
    splicebed = dd.get_junction_bed(data)
    if not file_exists(splicebed):
        logger.info("Splice junction BED file not found, skipping filtering of "
                    "variants closed to splice junctions.")
        return vrn_file
    spliceslop = get_padded_bed_file(vrn_dir, splicebed, SJ_BP_MASK, data)
    out_file = os.path.splitext(vrn_file)[0] + "-junctionfiltered.vcf.gz"
    if file_exists(out_file):
        return out_file
    with file_transaction(data, out_file) as tx_out_file:
        out_base = os.path.splitext(tx_out_file)[0]
        logger.info("Removing variants within %d bases of splice junctions listed in %s from %s. " % (SJ_BP_MASK, spliceslop, vrn_file))
        pybedtools.BedTool(vrn_file).intersect(spliceslop, wa=True, header=True, v=True).saveas(out_base)
        tx_out_file = vcfutils.bgzip_and_index(out_base, dd.get_config(data))
    return out_file
