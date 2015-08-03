"""Identify windows with very high depth for potential filtering.

In non-targeted experiments, high depth regions are often due to collapsed repeats
or other structure which can create long run times and incorrect results in
small and structural variant calling.
"""
import os
import shutil
import subprocess
import sys

import numpy
import yaml

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do

def _get_files(data):
    work_bam = dd.get_align_bam(data) or dd.get_work_bam(data)
    out_file = "%s-highdepth.bed" % utils.splitext_plus(work_bam)[0]
    stats_file = "%s-stats.yaml" % utils.splitext_plus(out_file)[0]
    return work_bam, out_file, stats_file

def identify(data):
    """Identify high depth regions in the alignment file for potential filtering.
    """
    high_multiplier = 20
    sample_size = int(1e6)
    high_percentage = 25.0
    min_coverage = 10
    window_size = 250
    work_bam, out_file, stats_file = _get_files(data)
    if not os.path.exists(out_file):
        cores = dd.get_num_cores(data)
        with file_transaction(data, out_file) as tx_out_file:
            tx_raw_file = "%s-raw%s" % utils.splitext_plus(tx_out_file)
            py_cl = os.path.join(os.path.dirname(sys.executable), "py")
            cmd = ("sambamba depth window -t {cores} -c {min_coverage} "
                   "--window-size {window_size} {work_bam} "
                   "| head -n {sample_size} "
                   """| cut -f 5 | {py_cl} -l 'numpy.median([float(x) for x in l if not x.startswith("mean")])'""")
            median_cov = float(subprocess.check_output(cmd.format(**locals()), shell=True))
            if not numpy.isnan(median_cov):
                high_thresh = int(high_multiplier * median_cov)
                cmd = ("sambamba depth window -t {cores} -c {median_cov} "
                       "--window-size {window_size} -T {high_thresh} {work_bam} "
                       "| {py_cl} -fx 'float(x.split()[5]) >= {high_percentage} "
                       """if not x.startswith("#") else None' """
                       "| cut -f 1-3,7 > {tx_raw_file} ")
                do.run(cmd.format(**locals()), "Identify high coverage regions")
                with open(stats_file, "w") as out_handle:
                    yaml.safe_dump({"median_cov": median_cov}, out_handle,
                                   allow_unicode=False, default_flow_style=False)
            else:
                with open(tx_raw_file, "w") as out_handle:
                    out_handle.write("")
            if utils.file_exists(tx_raw_file):
                cmd = "bedtools merge -i {tx_raw_file} -c 4 -o distinct > {tx_out_file}"
                do.run(cmd.format(**locals()), "Clean up raw coverage file")
            else:
                shutil.move(tx_raw_file, tx_out_file)
    return out_file if os.path.exists(out_file) else None

def get_median_coverage(data):
    stats_file = _get_files(data)[-1]
    if not utils.file_exists(stats_file):
        return 0
    else:
        with open(stats_file) as in_handle:
            stats = yaml.safe_load(in_handle)
        return stats["median_cov"]
