"""
Some reports for the analysis
"""
import shutil
import os
import sys
import copy
from contextlib import closing
from collections import Counter
import pandas as pd
import math

import pysam
import scipy.stats as stats

from bcbio.utils import splitext_plus, file_exists, safe_makedir, chdir, append_stem
from bcbio.log import logger
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio import bam
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.pipeline.alignment import get_aligner_index


def _run_meth_extractor(bam_in, sample, workdir, index_dir, config):
    """Run bismark_methylation_extractor command"""
    bismark = config_utils.get_program("bismark_methylation_extractor", config)
    cores = config["algorithm"].get('cores', 1)
    memory = config["algorithm"].get('mem', 5)
    bam_in = bam.sort(bam_in, config, order="queryname")
    cmd = "{bismark} --no_overlap --comprehensive --cytosine_report --genome_folder {index_dir} --merge_non_CpG --multicore {cores} --buffer_size {memory}G --bedGraph --gzip {bam_in}"
    out_dir = os.path.join(workdir, sample)
    mbias_file = os.path.join(out_dir, os.path.basename(splitext_plus(bam_in)[0]) + '.M-bias.txt')
    if not file_exists(mbias_file):
        with tx_tmpdir() as tx_dir:
            with chdir(tx_dir):
                do.run(cmd.format(**locals()), "bismark_methylation_extractor  in %s" % bam_in)
                shutil.move(tx_dir, out_dir)
    assert os.path.exists(mbias_file), "mbias report doesn't exists:%s" % mbias_file
    return mbias_file


def _run_report(bam_in, bam_report, sample, biasm_file, workdir, config):
    """
    Run bismark2report command
    """
    bismark = config_utils.get_program("bismark2report", config)
    cmd = "{bismark} --alignment_report {bam_report} -o {tx_out} --mbias_report {biasm_file}"
    out_dir = os.path.join(workdir, sample)
    out_file = os.path.join(out_dir, sample + '.html')
    with chdir(out_dir):
        if not file_exists(out_file):
            with file_transaction(out_file) as tx_out:
                do.run(cmd.format(**locals()), "bismarkr2report  in %s" % bam_in)
    return out_dir

def _bismark_calling(data):
    workdir = safe_makedir(os.path.join(dd.get_work_dir(data), "cpg"))
    config = data["config"]
    sample = dd.get_sample_name(data)
    index_dir = get_aligner_index('bismark', data)
    biasm_file = _run_meth_extractor(data["work_bam"], sample, workdir, index_dir, config)
    data['bismark_report'] = _run_report(data["work_bam"], data["bam_report"], sample, biasm_file, workdir, config)
    return data

def _bsmap_calling(data):
    sample = dd.get_sample_name(data)
    workdir = safe_makedir(os.path.join(dd.get_work_dir(data), "cpg_split", sample))
    config = data["config"]
    ref = dd.get_sam_ref(data)
    work_bam = dd.get_work_bam(data)
    python = os.path.join(os.path.dirname(sys.executable), "python")
    methratio = config_utils.get_program("methratio.py", config)
    cmd = ("{python} {methratio} -g -n -u -p -r -m 5 --chr={chrom} --ref={ref} {work_bam} >> {out_tx}")
    chrom = data["chr_to_run"]

    out_file = os.path.join(workdir, "methyratios_%s.txt" % chrom)
    if not file_exists(out_file):
        with file_transaction(out_file) as out_tx:
            do.run(cmd.format(**locals()), "Extract methylation for: %s" % sample)
    data["cpg_file"] = out_file
    return data

def calling(data):
    if dd.get_aligner(data) == "bismark":
        data = _bismark_calling(data)
    if dd.get_aligner(data) == "bsmap":
        data = _bsmap_calling(data)
    return [[data]]

# All these functions were very specific to a consult. Probably remove is
# the best option since they are not run for the general pipeline.
def parallel_calling(data, run_parallel):
    """This is needed only if running methylated veruss hidroxy-methulated"""
    out = []
    for sample in data:
        work_bam = dd.get_work_bam(sample[0])
        with closing(pysam.Samfile(work_bam, "rb")) as pysam_work_bam:
            chroms = pysam_work_bam.references
            for chrom in chroms:
                new_sample = copy.deepcopy(sample)
                if chrom.find("_") > -1:
                    continue
                new_sample[0]['chr_to_run'] = chrom
                out.append(new_sample)
    out = run_parallel("cpg_calling", out)
    for sample in out:
        phenotype = dd.get_phenotype(sample[0])
        batch = dd.get_batch(sample[0])
        if phenotype == "mC":
            for sample2 in out:
                if batch in dd.get_batch(sample2[0]) and dd.get_phenotype(sample2[0]) == "hmC":
                    if sample[0]["chr_to_run"] == sample2[0]["chr_to_run"]:
                        sample[0]["control"] = sample2[0]["cpg_file"]
                        break
    out = run_parallel("cpg_processing", out)
    for sample in data:
        sample[0]["cpg_split"] = []
        sample[0]["hmc_split"] = []
        name = dd.get_sample_name(sample[0])
        for chunck in out:
            if name == dd.get_sample_name(chunck[0]):
                sample[0]["cpg_split"].append(chunck[0]["cpg_file"])
                if "hmc_file" in chunck[0]:
                    sample[0]["hmc_split"].append(chunck[0]["hmc_file"])
    #  run_parallel("cpg_stats", data)

def _sync_pos(handle, tag):
    for line in handle:
        cols = line.strip().split("\t")
        info = cols[4:]
        if cols[3] != "CG":
            continue
        pos = int(cols[1])
        ratio = [int(float(cols[5])), int(cols[6])]
        if tag <= pos:
            return [handle, {"pos": pos, "counts": ratio, "info": info, "ratio": float(cols[4])}]
    return [None, None]

def _call_hmc(mc, hmc):
    n1, n2 = mc[0] - mc[1], mc[1]
    n1 = max(0, n1)
    m1, m2 = hmc[0] - hmc[1], hmc[1]
    m1 = max(0, m1)
    oddsratio, pvalue = stats.fisher_exact([[n1, n2], [m1, m2]])
    return pvalue

def cpg_postprocessing(data):
    mC = data["cpg_file"]
    if not "control" in data:
        return [[data]]
    hmC = data["control"]
    out_file = append_stem(mC, "_hmC")
    pos = 0
    pos_hmC = 0
    data["hmc_file"] = out_file
    if file_exists(out_file):
        return [[data]]
    logger.debug("processing %s versus %s" % (mC, hmC))
    with file_transaction(out_file) as out_tx:
        with open(out_tx, "w") as out_handle:
            with open(mC) as mC_h:
                with open(hmC) as hmC_h:
                    for line in mC_h:
                        cols = line.strip().split("\t")
                        if cols[3] != "CG":
                            continue
                        pos = int(cols[1])
                        counts = [int(float(cols[5])), int(cols[6])]
                        if pos < pos_hmC:
                            continue
                        elif pos > pos_hmC:
                            hmC_h, hmC = _sync_pos(hmC_h, pos)
                            if not hmC_h:
                                break
                            pos_hmC = hmC["pos"]
                        if counts[0] < 9 or hmC["counts"][0] < 9:
                            continue
                        if pos == hmC["pos"]:
                            pvalue = _call_hmc(counts, hmC["counts"])
                            print >>out_handle, "%s\t%s\t%s" % (line.strip(), "\t".join(hmC["info"]), pvalue)
    return [[data]]

def hmc_stats(sample):
    dt = Counter()
    work_dir = dd.get_work_dir(sample)
    sample_name = dd.get_sample_name(sample)
    out = os.path.join(work_dir, "cpg_split", sample_name, "%s_hmC.tsv" % sample_name)
    cmd = "cat %s > {out}" % " ".join(sample["hmc_split"])
    do.run(cmd.format(**locals()), "")

def cpg_stats(sample):
    dtdepth = Counter()
    dtratio = Counter()
    work_dir = dd.get_work_dir(sample)
    sample_name = dd.get_sample_name(sample)
    depth_out = os.path.join(work_dir, "cpg_split", sample_name, "depth.tsv")
    ratio_out = os.path.join(work_dir, "cpg_split", sample_name, "ratio.tsv")
    hmc_out = os.path.join(work_dir, "cpg_split", sample_name, "%s_hmc.tsv.gz" % sample_name)
    hmc_files = " ".join(sample["hmc_split"])
    with file_transaction(hmc_out) as tx_out:
        header = " ".join(["chr", "pos", "strand", "context", "ratio", "eff_CT_counts",
                           "C_counts", "CT_counts", "rev_G_counts", "rev_GA_counts",
                           "CI_lover", "CI_upper", "ox_ratio", "ox_eff_CT_counts",
                           "ox_C_counts", "ox_CT_counts", "ox_rev_G_counts",
                           "ox_rev_GA_counts", "ox_CI_lower", "ox_CI_upper", "pvalue"])
        cmd = "cat <(echo {header} | sed 's/ /\t/g') {hmc_files} | gzip -c > {tx_out}"
        if not file_exists(hmc_out):
            do.run(cmd.format(**locals()), "Merging %s" % sample_name)
    work_bam = dd.get_work_bam(sample)
    if not file_exists(depth_out):
        for cpg_file in sample["cpg_split"]:
            logger.debug("Reading %s of sample %s" % (cpg_file, sample_name))
            if file_exists(cpg_file):
                 with open(cpg_file) as in_handle:
                     for line in in_handle:
                         cols = line.strip().split("\t")
                         if cols[3] == "CG":
                             ratio = int(float(cols[4]) * 100)
                             dtratio[ratio] += 1
                             depth = int(math.ceil(float(cols[5]))) if float(cols[5]) < 50 else 50
                             dtdepth[depth] += 1
        pd.DataFrame(dtdepth, index=[1]).to_csv(depth_out, sep="\t")
        pd.DataFrame(dtratio, index=[1]).to_csv(ratio_out, sep="\t")
    # calculate mlml
    if not hmc_files:
        return None
    out_dir = safe_makedir(os.path.join(work_dir, "mlml", sample_name))
    mlml_out = os.path.join(out_dir, "%s_mlml.txt.gz" % sample_name)
    if not file_exists(mlml_out) and file_exists(hmc_out):
        with chdir(out_dir):
            with file_transaction(mlml_out) as tx_out:
                tx_out_1 = "%s_noheader" % tx_out
                tx_out_2 = "%s_alone" % tx_out
                cmd = " ".join(["zcat %s | sed -e '1d' | awk " % hmc_out, ''' '{rounded = sprintf("%d", $14);print $1"\\t"$2"\t"$3"\\tCpG\\t"$13"\\t"rounded}' ''', "> %s_ox.txt" % sample_name])
                do.run(cmd, "Creating OX input for %s" % sample_name)
                cmd = " ".join(["zcat %s | sed -e '1d' | awk " % hmc_out, ''' '{rounded = sprintf("%d", $6);print $1"\\t"$2"\\t"$3"\\tCpG\\t"$5"\\t"rounded}' ''', "> %s_bs.txt" % sample_name])
                do.run(cmd, "Creating BS input for %s" % sample_name)
                cmd = ("mlml -o {tx_out_1} -u {sample_name}_bs.txt -m {sample_name}_ox.txt -v").format(**locals())
                do.run(cmd, "Run MLML with %s" % sample_name)
                cmd = ("cat <(echo chrom start end mC hmC C conflicts | sed 's/ /\t/g') {tx_out_1} | gzip -c > {tx_out_2} ").format(**locals())
                do.run(cmd, "")
                tx_out_1 = "%s.woFDR.gz" % tx_out
                cmd = ("paste <(zcat {hmc_out}) <(zcat {tx_out_2}) | gzip -c > {tx_out}").format(**locals())
                do.run(cmd, "Merge data for %s" % sample_name)

    merge_out = [os.path.join(out_dir, "%s_merged.txt.gz" % sample_name), os.path.join(out_dir, "%s_merged_pass.txt.gz" % sample_name)]
    if not file_exists(merge_out[0]):
        with file_transaction(merge_out) as tx_outs:
            tx_out, tx_out_pass = tx_outs
            df = pd.read_csv(mlml_out, sep="\t")
            import statsmodels.sandbox.stats.multicomp
            df["fdr"] = statsmodels.sandbox.stats.multicomp.fdrcorrection0(df["pvalue"])[1]
            df_p_pass = df[df.fdr<0.05]
            logger.debug("Pass FDR 5 pct in %s:%s " % (sample_name, float(df_p_pass.shape[1])/float(df.shape[1])))
            df.to_csv(tx_out, sep="\t")
            df_p_pass.to_csv(tx_out_pass, sep="\t")
