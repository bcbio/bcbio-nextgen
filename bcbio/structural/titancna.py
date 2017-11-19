"""TitanCNA: Subclonal CNV calling and loss of heterogeneity in cancer.

https://github.com/gavinha/TitanCNA
"""
import csv
import os
import shutil

import pandas as pd

from bcbio import utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import vcfutils

def run(items):
    from bcbio import heterogeneity
    paired = vcfutils.get_paired(items)
    if not paired:
        logger.info("Skipping TitanCNA; no somatic tumor calls in batch: %s" %
                    " ".join([dd.get_sample_name(d) for d in items]))
        return items
    work_dir = _sv_workdir(paired.tumor_data)
    cn_file = _titan_cn_file(dd.get_normalized_depth(paired.tumor_data), work_dir, paired.tumor_data)
    het_file = _titan_het_file(heterogeneity.get_variants(paired.tumor_data), work_dir, paired)
    if _should_run(het_file):
        ploidy_outdirs = []
        for ploidy in [2, 3, 4]:
            for num_clusters in [1, 2, 3]:
                out_dir = _run_titancna(cn_file, het_file, ploidy, num_clusters, work_dir, paired.tumor_data)
                ploidy_outdirs.append((ploidy, out_dir))
        solution_file = _run_select_solution(ploidy_outdirs, work_dir, paired.data)
        print(solution_file)
    else:
        logger.info("Skipping TitanCNA; not enough input data: %s" %
                    " ".join([dd.get_sample_name(d) for d in items]))
    return items

def _should_run(het_file):
    """Check for enough input data to proceed with analysis.
    """
    has_hets = False
    with open(het_file) as in_handle:
        for i, line in enumerate(in_handle):
            if i > 1:
                has_hets = True
                break
    return has_hets

def _run_select_solution(ploidy_outdirs, work_dir, data):
    """Select optimal
    """
    out_file = os.path.join(work_dir, "optimalClusters.txt")
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            ploidy_inputs = " ".join(["--ploidyRun%s=%s" % (p, d) for p, d in ploidy_outdirs])
            cmd = "titanCNA_selectSolution.R {ploidy_inputs} --outFile={tx_out_file}"
            do.run(cmd.format(**locals()), "TitanCNA: select optimal solution")
    return out_file

def _run_titancna(cn_file, het_file, ploidy, num_clusters, work_dir, data):
    """Run titanCNA wrapper script on given ploidy and clusters.
    """
    sample = dd.get_sample_name(data)
    cores = dd.get_num_cores(data)
    export_cmd = utils.get_R_exports()
    ploidy_dir = utils.safe_makedir(os.path.join(work_dir, "run_ploidy%s" % ploidy))
    cluster_dir = "%s_cluster%02d" % (sample, num_clusters)
    out_dir = os.path.join(ploidy_dir, cluster_dir)
    if not utils.file_exists(out_dir):
        with tx_tmpdir(data) as tmp_dir:
            cmd = ("{export_cmd} && titanCNA.R --id {sample} --hetFile {het_file} --cnFile {cn_file} "
                   "--numClusters {num_clusters} --ploidy {ploidy} --numCores {cores} --outDir {tmp_dir}")
            do.run(cmd.format(**locals()), "TitanCNA CNV detection: ploidy %s, cluster %s" % (ploidy, num_clusters))
            shutil.move(os.path.join(tmp_dir, cluster_dir), out_dir)
    return out_dir

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(dd.get_work_dir(data), "structural",
                                           dd.get_sample_name(data), "titancna"))

def _titan_het_file(vrn_files, work_dir, paired):
    assert vrn_files, "Did not find compatible variant calling files for TitanCNA inputs"
    from bcbio.heterogeneity import bubbletree
    class OutWriter:
        def __init__(self, out_handle):
            self.writer = csv.writer(out_handle, dialect="excel-tab")
        def write_header(self):
            self.writer.writerow(["Chr", "Position", "Ref", "RefCount", "Nref", "NrefCount", "NormQuality"])
        def write_row(self, rec, stats):
            self.writer.writerow([rec.chrom, rec.pos, rec.ref, stats["tumor"]["depth"] - stats["tumor"]["alt"],
                                  rec.alts[0], stats["tumor"]["alt"], rec.qual])
    return bubbletree.prep_vrn_file(vrn_files[0]["vrn_file"], vrn_files[0]["variantcaller"],
                                    work_dir, paired, OutWriter)

def _titan_cn_file(cnr_file, work_dir, data):
    """Convert generallized CNVkit normalized input into TitanCNA ready format.
    """
    out_file = os.path.join(work_dir, "%s.cn" % (utils.splitext_plus(os.path.basename(cnr_file))[0]))

    if not utils.file_uptodate(out_file, cnr_file):
        with file_transaction(data, out_file) as tx_out_file:
            iterator = pd.read_table(cnr_file, iterator=True)
            with open(tx_out_file, "w") as handle:
                for chunk in iterator:
                    chunk = chunk[["chromosome", "start", "end", "log2"]]
                    chunk.columns = ["chrom", "start", "end", "logR"]
                    chunk['start'] += 1
                    chunk.to_csv(handle, mode="a", sep="\t", index=False)
    return out_file
