"""PureCN: Copy number calling and SNV classification using targeted short read sequencing

https://github.com/lima1/PureCN
"""
import os
import re
import shutil
import subprocess

import pandas as pd
import toolz as tz

from bcbio import utils
from bcbio.heterogeneity import chromhacks
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.distributed.transaction import file_transaction
from bcbio.heterogeneity import loh
from bcbio.provenance import do
from bcbio.variation import germline, vcfutils
from bcbio.structural import cnvkit, gatkcnv, regions

def run(items):
    paired = vcfutils.get_paired(items)
    # paired is PairedInfo of one T/N pair (or just T) - named tuple, paired.tumor_config
    if not paired:
        logger.info("Skipping PureCN; no somatic tumor calls in batch: %s" %
                    " ".join([dd.get_sample_name(d) for d in items]))
        return items
    work_dir = _sv_workdir(paired.tumor_data)
    normaldb = tz.get_in(["algorithm", "background", "cnv_reference", "purecn_normaldb"], paired.tumor_config)
    # the right way of running purecn is with normaldb
    if normaldb:
        purecn_out = _run_purecn_normaldb(paired, work_dir)
        purecn_out = _run_purecn_dx(purecn_out, paired)
    else:
        purecn_out = _run_purecn(paired, work_dir)
    out = []
    if paired.normal_data:
        out.append(paired.normal_data)
    if purecn_out:
        purecn_out["variantcaller"] = "purecn"
        if "loh" in purecn_out:
            from bcbio.structural import titancna
            purecn_out["vrn_file"] = titancna.to_vcf(purecn_out["loh"], "PureCN", _get_header, _loh_to_vcf,
                                                     paired.tumor_data, sep=",")
            purecn_out["lohsummary"] = loh.summary_status(purecn_out, paired.tumor_data)
        if "sv" not in paired.tumor_data:
            paired.tumor_data["sv"] = []
        paired.tumor_data["sv"].append(purecn_out)
    out.append(paired.tumor_data)
    return out

def _run_purecn_normaldb(paired, out):
    """Run PureCN with normaldb and native segmentation
       paired is one t/n pair or only """
    sample = utils.to_single_data(paired.tumor_data)
    bed_file = tz.get_in(["config", "algorithm", "purecn_bed_ready"], sample)
    sample_name = dd.get_sample_name(sample)
    work_dir = _sv_workdir(sample)
    rscript = utils.Rscript_cmd()
    purecn_r = utils.R_package_script("PureCN", "extdata/PureCN.R")
    intervals = tz.get_in(["config", "algorithm", "purecn_bed_ready"], sample)
    bam_file = dd.get_align_bam(sample)
    # termline and somatic - just annotated and filters assigned
    variants_vcf =  tz.get_in(["variants"], sample)[0].get("germline")
    # in a T/N case, there is no germline file - vrn file with all variants
    if not variants_vcf:
        variants_vcf = tz.get_in(["variants"], sample)[0].get("vrn_file")
    normaldb = tz.get_in(["config", "algorithm", "background", "cnv_reference", "purecn_normaldb"], sample)
    mappingbiasfile = tz.get_in(["config", "algorithm", "background", "cnv_reference", "purecn_mapping_bias"], sample)
    sample_coverage = tz.get_in(["depth", "bins", "purecn"], sample)
    simple_repeat_bed = dd.get_variation_resources(sample)["simple_repeat"]
    result_file = os.path.join(work_dir, sample_name + ".rds")
    genome = dd.get_genome_build(sample)
    cmd = [ rscript, purecn_r,
            "--out", work_dir,
            "--tumor", sample_coverage,
            "--sampleid", sample_name,
            "--vcf", variants_vcf,
            "--normaldb", normaldb,
            "--mappingbiasfile", mappingbiasfile,
            "--intervals", intervals,
            "--snpblacklist", simple_repeat_bed,
            "--genome", genome,
            "--force",
            "--postoptimize",
            "--seed", "123",
            "--bootstrapn", "500",
            "--cores", dd.get_num_cores(sample)]
    resources = config_utils.get_resources("purecn", sample)
    if "options" in resources:
        cmd += [str(x) for x in resources.get("options", [])]
    # it is not recommended to use matched normal sample in PureCN analysis,
    # because then it skips PON coverage normalization and denoising steps!
    # but still, if it is supplied, we useit
    if paired.normal_data:
        normal_sample = utils.to_single_data(paired.normal_data)
        if normal_sample:
            normal_coverage = tz.get_in(["depth", "bins", "purecn"], normal_sample)
            cmd.extend(["--normal", normal_coverage])
    if not os.path.exists(result_file):
        try:
            cmd_line = "export R_LIBS_USER=%s && %s && %s" % (utils.R_sitelib(env = "base"),
                                                              utils.get_R_exports(env = "base"),
                                                              " ".join([str(x) for x in cmd]))
            do.run(cmd_line, "PureCN copy number calling")
            logger.debug("Saved PureCN output to " + work_dir)
        except subprocess.CalledProcessError as msg:
            logger.info("PureCN failed")
    out_base, out, all_files  = _get_purecn_files(paired, work_dir, require_exist = True)
    return out

def _run_purecn_dx(out, paired):
    """Extract signatures and mutational burdens from PureCN rds file."""
    out_base, out, all_files = _get_purecn_dx_files(paired, out, require_exist = True)
    rscript = utils.Rscript_cmd()
    purecndx_r = utils.R_package_script("PureCN", "extdata/Dx.R")
    simple_repeat_bed = dd.get_variation_resources(paired.tumor_data)["simple_repeat"]
    callable_bed = dd.get_sample_callable(paired.tumor_data)
    if not utils.file_uptodate(out["mutation_burden"], out["rds"]):
        with file_transaction(paired.tumor_data, out_base) as tx_out_base:
            cmd = [rscript, purecndx_r, 
                   "--rds", out["rds"], 
                   "--callable", callable_bed,
                   "--signatures",
                   "--exclude", simple_repeat_bed,
                   "--out", tx_out_base]
            do.run(cmd, "PureCN Dx mutational burden and signatures")
            for f in all_files:
                if os.path.exists(os.path.join(os.path.dirname(tx_out_base), f)):
                    shutil.move(os.path.join(os.path.dirname(tx_out_base), f),
                                os.path.join(os.path.dirname(out_base), f))
    return out

def _get_purecn_dx_files(paired, out, require_exist = False):
    """Retrieve files generated by PureCN_Dx"""
    out_base = utils.splitext_plus(out["rds"])[0]
    all_files = []
    for key, ext in [[("mutation_burden",), "_mutation_burden.csv"],
                     [("plot", "signatures"), "_signatures.pdf"],
                     [("signatures",), "_signatures.csv"],
                     [("chrom_instability",), "_cin.csv"]]:
        cur_file = f"{out_base}{ext}"
        if not require_exist or os.path.exists(cur_file):
            out = tz.update_in(out, key, lambda x: cur_file)
            all_files.append(os.path.basename(cur_file))
    return out_base, out, all_files

def _run_purecn(paired, work_dir):
    """Run PureCN.R wrapper with pre-segmented CNVkit or GATK4 inputs.
    """
    segfns = {"cnvkit": _segment_normalized_cnvkit, "gatk-cnv": _segment_normalized_gatk}
    out_base, out, all_files = _get_purecn_files(paired, work_dir)
    failed_file = out_base + "-failed.log"
    cnr_file = tz.get_in(["depth", "bins", "normalized"], paired.tumor_data)
    if not utils.file_uptodate(out["rds"], cnr_file) and not utils.file_exists(failed_file):
        cnr_file, seg_file = segfns[cnvkit.bin_approach(paired.tumor_data)](cnr_file, work_dir, paired)
        from bcbio import heterogeneity
        vcf_file = heterogeneity.get_variants(paired.tumor_data, include_germline=False)[0]["vrn_file"]
        vcf_file = germline.filter_to_pass_and_reject(vcf_file, paired, out_dir=work_dir)
        with file_transaction(paired.tumor_data, out_base) as tx_out_base:
            # Use UCSC style naming for human builds to support BSgenome
            genome = ("hg19" if dd.get_genome_build(paired.tumor_data) in ["GRCh37", "hg19"]
                      else dd.get_genome_build(paired.tumor_data))
            rscript = utils.Rscript_cmd()
            purecn_r = utils.R_package_script("PureCN", "extdata/PureCN.R")
            cmd = [rscript, purecn_r, "--seed", "42", "--out", tx_out_base, 
                   "--rds", "%s.rds" % tx_out_base,
                   "--sampleid", dd.get_sample_name(paired.tumor_data),
                   "--genome", genome,
                   "--vcf", vcf_file, "--tumor", cnr_file,
                   "--segfile", seg_file, "--funsegmentation", "Hclust", "--maxnonclonal", "0.3"]
            if dd.get_num_cores(paired.tumor_data) > 1:
                cmd += ["--cores", str(dd.get_num_cores(paired.tumor_data))]
            try:
                cmd = "export R_LIBS_USER=%s && %s && %s" % (utils.R_sitelib(env="base"),
                                                             utils.get_R_exports(env="base"),
                                                             " ".join([str(x) for x in cmd]))
                do.run(cmd, "PureCN copy number calling")
            except subprocess.CalledProcessError as msg:
                if _allowed_errors(str(msg)):
                    logger.info("PureCN failed to find solution for %s: skipping" %
                                dd.get_sample_name(paired.tumor_data))
                    with open(failed_file, "w") as out_handle:
                        out_handle.write(str(msg))
                else:
                    logger.exception()
                    raise
            for f in all_files:
                if os.path.exists(os.path.join(os.path.dirname(tx_out_base), f)):
                    shutil.move(os.path.join(os.path.dirname(tx_out_base), f),
                                os.path.join(os.path.dirname(out_base), f))
    out = _get_purecn_files(paired, work_dir, require_exist=True)[1]
    return out if (out.get("rds") and os.path.exists(out["rds"])) else None

def _allowed_errors(msg):
    allowed = ["Could not find valid purity and ploidy solution.",
               "Cannot find valid purity/ploidy solution",
               "None of the variants in provided VCF passed filtering."]
    return any([len(re.findall(m, msg)) > 0 for m in allowed])

def _segment_normalized_gatk(cnr_file, work_dir, paired):
    """Segmentation of normalized inputs using GATK4, converting into standard input formats.
    """
    work_dir = utils.safe_makedir(os.path.join(work_dir, "gatk-cnv"))
    seg_file = gatkcnv.model_segments(cnr_file, work_dir, paired)["seg"]
    std_seg_file = seg_file.replace(".cr.seg", ".seg")
    if not utils.file_uptodate(std_seg_file, seg_file):
        with file_transaction(std_seg_file) as tx_out_file:
            df = pd.read_csv(seg_file, sep="\t", comment="@", header=0,
                             names=["chrom", "loc.start", "loc.end", "num.mark", "seg.mean"])
            df.insert(0, "ID", [dd.get_sample_name(paired.tumor_data)] * len(df))
            df.to_csv(tx_out_file, sep="\t", header=True, index=False)
    std_cnr_file = os.path.join(work_dir, "%s.cnr" % dd.get_sample_name(paired.tumor_data))
    if not utils.file_uptodate(std_cnr_file, cnr_file):
        with file_transaction(std_cnr_file) as tx_out_file:
            logdf = pd.read_csv(cnr_file, sep="\t", comment="@", header=0,
                                names=["chrom", "start", "end", "log2"])
            covdf = pd.read_csv(tz.get_in(["depth", "bins", "antitarget"], paired.tumor_data),
                                sep="\t", header=None,
                                names=["chrom", "start", "end", "orig.name", "depth", "gene"])
            df = pd.merge(logdf, covdf, on=["chrom", "start", "end"])
            del df["orig.name"]
            df = df[["chrom", "start", "end", "gene", "log2", "depth"]]
            df.insert(6, "weight", [1.0] * len(df))
            df.to_csv(tx_out_file, sep="\t", header=True, index=False)
    return std_cnr_file, std_seg_file

def _segment_normalized_cnvkit(cnr_file, work_dir, paired):
    """Segmentation of normalized inputs using CNVkit.
    """
    cnvkit_base = os.path.join(utils.safe_makedir(os.path.join(work_dir, "cnvkit")),
                                dd.get_sample_name(paired.tumor_data))
    cnr_file = chromhacks.bed_to_standardonly(cnr_file, paired.tumor_data, headers="chromosome",
                                                include_sex_chroms=True,
                                                out_dir=os.path.dirname(cnvkit_base))
    cnr_file = _remove_overlaps(cnr_file, os.path.dirname(cnvkit_base), paired.tumor_data)
    seg_file = cnvkit.segment_from_cnr(cnr_file, paired.tumor_data, cnvkit_base)
    return cnr_file, seg_file

def _remove_overlaps(in_file, out_dir, data):
    """Remove regions that overlap with next region, these result in issues with PureCN.
    """
    out_file = os.path.join(out_dir, "%s-nooverlaps%s" % utils.splitext_plus(os.path.basename(in_file)))
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(in_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    prev_line = None
                    for line in in_handle:
                        if prev_line:
                            pchrom, pstart, pend = prev_line.split("\t", 4)[:3]
                            cchrom, cstart, cend = line.split("\t", 4)[:3]
                            # Skip if chromosomes match and end overlaps start
                            if pchrom == cchrom and int(pend) > int(cstart):
                                pass
                            else:
                                out_handle.write(prev_line)
                        prev_line = line
                    out_handle.write(prev_line)
    return out_file

def _get_purecn_files(paired, work_dir, require_exist=False):
    """Retrieve organized structure of PureCN output files."""
    sample_name = dd.get_sample_name(paired.tumor_data)
    out_base = os.path.join(work_dir, sample_name)
    out = {"plot": {}}
    all_files = []
    for plot in ["chromosomes", "local_optima", "segmentation", "summary"]:
        if plot == "summary":
            cur_file = f"{out_base}.pdf"
        else:
            cur_file = f"{out_base}_{plot}.pdf"
        if not require_exist or os.path.exists(cur_file):
            out["plot"][plot] = cur_file
            all_files.append(os.path.basename(cur_file))
    for key, ext in [["hetsummary", ".csv"],
                     ["dnacopy", "_dnacopy.seg"], 
                     ["genes", "_genes.csv"],
                     ["log", ".log"], 
                     ["loh", "_loh.csv"], 
                     ["rds", ".rds"],
                     ["variants", "_variants.csv"]]:
        cur_file = f"{out_base}{ext}"
        if not require_exist or os.path.exists(cur_file):
            out[key] = cur_file
            all_files.append(os.path.basename(cur_file))
    sample_name = dd.get_sample_name(paired.tumor_data)
    more_files = [sample_name + item for item in ["_amplification_pvalues.csv", "_chromosomes.pdf",
                                   ".csv", "_dnacopy.seg", "_genes.csv", "_local_optima.pdf",
                                   ".log", "_loh.csv", ".pdf", ".rds", "_segmentation.pdf",
                                   "_variants.csv"]]
    for purecn_file in more_files:
        purecn_file_path = os.path.join(work_dir, purecn_file)
        if not require_exist or os.path.exists(purecn_file_path):
            out[purecn_file] = purecn_file_path
            all_files.append(purecn_file_path)
    return out_base, out, all_files

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(dd.get_work_dir(data), "structural",
                                           dd.get_sample_name(data), "purecn"))

# ## VCF output

def _get_header(in_handle):
    return in_handle.readline().strip().split(","), in_handle

def _loh_to_vcf(cur):
    """Convert LOH output into standardized VCF."""
    # PureCN 1.14 outputs segments without informative SNPs, skip those
    # see https://github.com/lima1/PureCN/blob/ed7d10c7ca578bc7d1aabef86893c23ddddf79dc/NEWS#L92-L94
    if cur["C"] == "NA" or cur["M"] == "NA":
        return None
    cn = int(float(cur["C"]))
    minor_cn = int(float(cur["M"]))
    if cur["type"].find("LOH") > -1:
        svtype = "LOH"
    elif cn > 2:
        svtype = "DUP"
    elif cn < 1:
        svtype = "DEL"
    else:
        svtype = None
    if svtype:
        # end could be 100.5
        start = int(float(cur["start"]))
        end = int(float(cur["end"]))
        info = [f"SVTYPE={svtype}", 
                f"END={end}",
                f"SVLEN={end-start+1}",
                f"CN={cn}",
                f"MajorCN={cn - minor_cn}",
                f"MinorCN={minor_cn}"]
        return [cur["chr"], cur["start"], ".", "N", "<%s>" % svtype, ".", ".",
                ";".join(info), "GT", "0/1"]

def process_intervals(data):
    """Prepare intervals file"""
    bed_file = regions.get_sv_bed(data)
    if not bed_file:
         bed_file = bedutils.clean_file(dd.get_variant_regions(data), data)
    if not bed_file:
        return None

    basename = os.path.splitext(bed_file)[0]
    ready_file = basename + ".txt"
    if os.path.exists(ready_file):
        return ready_file
    optimized_bed = basename + ".optimized.bed"
    rscript = utils.Rscript_cmd("base")
    interval_file_r = utils.R_package_script("base", "PureCN", "extdata/IntervalFile.R")
    ref_file = dd.get_ref_file(data)
    mappability_resource = dd.get_variation_resources(data)["purecn_mappability"]
    genome = dd.get_genome_build(data)
    cmd = [rscript, interval_file_r, "--infile", bed_file,
          "--fasta", ref_file,
          "--outfile", ready_file,
          "--offtarget",
          "--genome", genome,
          "--export", optimized_bed,
          "--mappability", mappability_resource]
    try:
        cmd_line = "export R_LIBS_USER=%s && %s && %s" % (utils.R_sitelib(env = "base"),
                                                     utils.get_R_exports(env = "base"),
                                                     " ".join([str(x) for x in cmd]))
        do.run(cmd_line, "PureCN intervals")
    except subprocess.CalledProcessError as msg:
        logger.info("PureCN failed to prepare intervals")
    logger.debug("Saved PureCN interval file into " + ready_file)
    return ready_file

def get_coverage(data):
    """Calculate coverage for a sample.bam, account for GC content
       data is single sample
    """
    data = utils.to_single_data(data)
    bed_file = tz.get_in(["config", "algorithm", "purecn_bed_ready"], data)
    sample_name = dd.get_sample_name(data)
    work_dir = _sv_workdir(data)
    rscript = utils.Rscript_cmd("base")
    coverage_r = utils.R_package_script("base", "PureCN", "extdata/Coverage.R")
    intervals = tz.get_in(["config", "algorithm", "purecn_bed_ready"], data)
    # PureCN resolves symlinks and the actual output PureCN coverage file name
    # is derived from the end bam not from bam_file
    bam_file = os.path.realpath(dd.get_align_bam(data))
    bam_name = os.path.basename(bam_file)
    (bname, ext) = os.path.splitext(bam_name)
    result_file = os.path.join(work_dir, bname + "_coverage_loess.txt.gz")
    if not os.path.exists(result_file):
        cmd = [rscript, coverage_r,
               "--outdir", work_dir,
               "--bam", bam_file,
               "--intervals", intervals]
        try:
            cmd_line = "export R_LIBS_USER=%s && %s && %s" % (utils.R_sitelib(env = "base"),
                                                              utils.get_R_exports(env = "base"),
                                                              " ".join([str(x) for x in cmd]))
            do.run(cmd_line, "PureCN coverage")
        except subprocess.CalledProcessError as msg:
            logger.info("PureCN failed to calculate coverage")
        logger.debug("Saved PureCN coverage files to " + result_file)
    return result_file

def create_normal_db(coverage_files_txt, snv_pon, out_dir, genome_build):
    """create normal db
       input: coverage files calculated by purecn for each sample
              snv_pon - mutect2 SNV PON
       output:
              mapping_bias_hg38.rds
              normalDB_hg38.rds
    """
    rscript = utils.Rscript_cmd("base")
    normaldb_r = utils.R_package_script("base", "PureCN", "extdata/NormalDB.R")
    cmd = [rscript, normaldb_r,
           "--outdir", out_dir,
           "--coveragefiles", coverage_files_txt,
           "--normal_panel" , snv_pon,
           "--genome", genome_build,
           "--force"]
    try:
        cmd_line = "export R_LIBS_USER=%s && %s && %s" % (utils.R_sitelib(env = "base"),
                                                          utils.get_R_exports(env = "base"),
                                                          " ".join([str(x) for x in cmd]))
        do.run(cmd_line, "PureCN normalDB")
    except subprocess.CalledProcessError as msg:
        logger.info("PureCN failed to create a normal db")

    return out_dir
