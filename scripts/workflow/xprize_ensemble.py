#!/usr/bin/env python
"""XPrize scoring workflow that converts input BAMs into consolidated Ensemble calls.

Automates the Ensemble approach described here (http://j.mp/VUbz9A) to prepare
a final set of reference haploid variant calls for X Prize scoring.

Usage:
  xprize_ensemble.py <sample name> <BAM file> <fosmid region BED file>
                     <base directory> <bcbio system YAML>
"""
import datetime
import os
import subprocess
import sys

import yaml

from bcbio import utils

def main(sample, bam_file, bed_file, base_dir, bcbio_config_file):
    configdir = utils.safe_makedir(os.path.join(base_dir, sample, "config"))
    out_config_file = os.path.join(configdir, "%s.yaml" % sample)
    callers = ["gatk", "freebayes", "samtools", "varscan"]
    out = {"fc_date": datetime.datetime.now().strftime("%y%m%d"),
           "fc_name": sample,
           "details": [{
               "files": [bam_file],
               "lane": 1,
               "description": sample,
               "analysis": "variant",
               "genome_build": "GRCh37",
               "algorithm": {
                   "aligner": False,
                   "recalibrate": False,
                   "realign": False,
                   "ploidy": 1,
                   "variantcaller": callers,
                   "quality_format": "Standard",
                   "variant_regions": bed_file,
                   "coverage_interval": "regional",
                   "ensemble": {
                       "format-filters": ["DP < 4"],
                       "classifiers": {
                           "balance": ["AD", "FS", "Entropy"],
                           "calling": ["ReadPosEndDist", "PL", "Entropy", "NBQ"]},
                       "classifier-params": {
                           "type": "svm"},
                       "trusted-pct": 0.65}}}]}
    with open(out_config_file, "w") as out_handle:
        yaml.dump(out, out_handle)

    workdir = utils.safe_makedir(os.path.join(base_dir, sample, "work"))
    os.chdir(workdir)
    subprocess.check_call(["bcbio_nextgen.py", bcbio_config_file, out_config_file, "-n", str(len(callers))])

if __name__ == "__main__":
    main(*sys.argv[1:])
