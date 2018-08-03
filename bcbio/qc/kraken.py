"""Contamination detection using Kraken.

https://ccb.jhu.edu/software/kraken/
"""
import os
import shutil
import toolz as tz

from bcbio import bam, install, utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.log import logger
from bcbio.provenance import do
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils

def run(_, data, out_dir):
    """Run kraken, generating report in specified directory and parsing metrics.
       Using only first paired reads.
    """
    # logger.info("Number of aligned reads < than 0.60 in %s: %s" % (dd.get_sample_name(data), ratio))
    logger.info("Running kraken to determine contaminant: %s" % dd.get_sample_name(data))
    # ratio = bam.get_aligned_reads(bam_file, data)
    out = out_stats = None
    db = tz.get_in(["config", "algorithm", "kraken"], data)
    if db and isinstance(db, (list, tuple)):
        db = db[0]
    kraken_cmd = config_utils.get_program("kraken", data["config"])
    if db == "minikraken":
        db = os.path.join(install._get_data_dir(), "genomes", "kraken", "minikraken")

    if not os.path.exists(db):
        logger.info("kraken: no database found %s, skipping" % db)
        return {"kraken_report": "null"}

    if not os.path.exists(os.path.join(out_dir, "kraken_out")):
        work_dir = os.path.dirname(out_dir)
        utils.safe_makedir(work_dir)
        num_cores = data["config"]["algorithm"].get("num_cores", 1)
        fn_file = data["files_orig"][0] if dd.get_save_diskspace(data) else data["files"][0]
        if fn_file.endswith("bam"):
            logger.info("kraken: need fastq files as input")
            return {"kraken_report": "null"}
        with tx_tmpdir(data) as tx_tmp_dir:
            with utils.chdir(tx_tmp_dir):
                out = os.path.join(tx_tmp_dir, "kraken_out")
                out_stats = os.path.join(tx_tmp_dir, "kraken_stats")
                cat = "zcat" if fn_file.endswith(".gz") else "cat"
                cl = ("{cat} {fn_file} | {kraken_cmd} --db {db} --quick "
                      "--preload --min-hits 2 "
                      "--threads {num_cores} "
                      "--output {out} --fastq-input /dev/stdin  2> {out_stats}").format(**locals())
                do.run(cl, "kraken: %s" % dd.get_sample_name(data))
                if os.path.exists(out_dir):
                    shutil.rmtree(out_dir)
                shutil.move(tx_tmp_dir, out_dir)
    metrics = _parse_kraken_output(out_dir, db, data)
    return metrics

def _parse_kraken_output(out_dir, db, data):
    """Parse kraken stat info comming from stderr,
       generating report with kraken-report
    """
    in_file = os.path.join(out_dir, "kraken_out")
    stat_file = os.path.join(out_dir, "kraken_stats")
    out_file = os.path.join(out_dir, "kraken_summary")
    kraken_cmd = config_utils.get_program("kraken-report", data["config"])
    classify = unclassify = None
    with open(stat_file, 'r') as handle:
        for line in handle:
            if line.find(" classified") > -1:
                classify = line[line.find("(") + 1:line.find(")")]
            if line.find(" unclassified") > -1:
                unclassify = line[line.find("(") + 1:line.find(")")]
    if os.path.getsize(in_file) > 0 and not os.path.exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cl = ("{kraken_cmd} --db {db} {in_file} > {tx_out_file}").format(**locals())
            do.run(cl, "kraken report: %s" % dd.get_sample_name(data))
    kraken = {"kraken_clas": classify, "kraken_unclas": unclassify}
    kraken_sum = _summarize_kraken(out_file)
    kraken.update(kraken_sum)
    return kraken

def _summarize_kraken(fn):
    """get the value at species level"""
    kraken = {}
    list_sp, list_value = [], []
    with open(fn) as handle:
        for line in handle:
            cols = line.strip().split("\t")
            sp = cols[5].strip()
            if len(sp.split(" ")) > 1 and not sp.startswith("cellular"):
                list_sp.append(sp)
                list_value.append(cols[0])
    kraken = {"kraken_sp": list_sp, "kraken_value": list_value}
    return kraken
