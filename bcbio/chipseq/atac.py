# add ["atac_metrics": {}] to the datadict

import os
import toolz as tz

from bcbio import utils
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.log import logger
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio import bam

def calculate_complexity_metrics(work_bam, data):
    """
    the work_bam should have duplicates marked but not removed
    mitochondrial reads should be removed 
    """
    bedtools = config_utils.get_program("bedtools", dd.get_config(data))
    work_dir = dd.get_work_dir(data)
    metrics_dir = os.path.join(work_dir, "metrics", "atac")
    utils.safe_makedir(metrics_dir)
    metrics_file = os.path.join(metrics_dir,
                                f"{dd.get_sample_name(data)}-atac-metrics.csv")
    if utils.file_exists(metrics_file):
        data = tz.assoc_in(data, ['atac', 'complexity_metrics_file'], metrics_file)
        return data
    # BAM file must be sorted by read name
    work_bam = bam.sort(work_bam, dd.get_config(data), order="queryname")
    with file_transaction(metrics_file) as tx_metrics_file:
        with open(tx_metrics_file, "w") as out_handle:
            out_handle.write("mt,m0,m1,m2\n")
        cmd = (f"{bedtools} bamtobed -bedpe -i {work_bam} | "
               "awk 'BEGIN{OFS=\"\\t\"}{print $1,$2,$4,$6,$9,$10}' | "
               "sort | "
               "uniq -c | "
               "awk 'BEGIN{mt=0;m0=0;m1=0;m2=0}($1==1){m1=m1+1} "
               "($1==2){m2=m2+1}{m0=m0+1}{mt=mt+$1}END{printf \"%d,%d,%d,%d\\n\", mt,m0,m1,m2}' >> "
               f"{tx_metrics_file}")
        message = f"Calculating ATAC-seq complexity metrics on {work_bam}, saving as {metrics_file}."
        do.run(cmd, message)
    data = tz.assoc_in(data, ['atac', 'complexity_metrics_file'], metrics_file)
    return data

def calculate_encode_complexity_metrics(data):
    metrics_file = tz.get_in(['atac', 'complexity_metrics_file'], data, None)
    if not metrics_file:
        return {}
    else:
        with open(metrics_file) as in_handle:
            header = next(in_handle).strip().split(",")
            values = next(in_handle).strip().split(",")
    raw_metrics = {h: int(v) for h, v in zip(header, values)}
    metrics = {"PBC1": raw_metrics["m1"] / raw_metrics["m0"],
               "NRF": raw_metrics["m0"] / raw_metrics["mt"]}
    if raw_metrics["m2"] == 0:
        PBC2 = 0
    else:
        PBC2 = raw_metrics["m1"] / raw_metrics["m2"]
    metrics["PBC2"] = PBC2
    return(metrics)
