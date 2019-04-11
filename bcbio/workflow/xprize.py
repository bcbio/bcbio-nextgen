"""XPrize scoring workflow that converts input BAMs into consolidated Ensemble calls.

Automates the Ensemble approach described here (http://j.mp/VUbz9A) to prepare
a final set of reference haploid variant calls for X Prize scoring.
"""
import argparse
import datetime
import os
import sys

import yaml

from bcbio import utils

class HelpArgParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args(args):
    parser = HelpArgParser(
        description="Automate ensemble variant approach for X Prize preparation")
    parser.add_argument("sample", help="Sample name")
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("bed_file", help="BED file of fosmid regions")
    parser.add_argument("base_dir", help="Base directory to process in")
    parser.add_argument("bcbio_config_file", help="bcbio system YAML config")
    args = parser.parse_args(args)
    return args

def get_fc_date(out_config_file):
    """Retrieve flowcell date, reusing older dates if refreshing a present workflow.
    """
    if os.path.exists(out_config_file):
        with open(out_config_file) as in_handle:
            old_config = yaml.safe_load(in_handle)
            fc_date = old_config["fc_date"]
    else:
        fc_date = datetime.datetime.now().strftime("%y%m%d")
    return fc_date

def setup(args):
    final_dir = utils.safe_makedir(os.path.join(args.base_dir, "ready"))
    configdir = utils.safe_makedir(os.path.join(args.base_dir, args.sample, "config"))
    out_config_file = os.path.join(configdir, "%s.yaml" % args.sample)
    callers = ["gatk", "freebayes", "samtools", "varscan"]
    out = {"fc_date": get_fc_date(out_config_file),
           "fc_name": args.sample,
           "upload": {"dir": final_dir},
           "details": [{
               "files": [args.bam_file],
               "lane": 1,
               "description": args.sample,
               "analysis": "variant",
               "genome_build": "GRCh37",
               "algorithm": {
                   "aligner": False,
                   "recalibrate": False,
                   "realign": False,
                   "ploidy": 1,
                   "variantcaller": callers,
                   "quality_format": "Standard",
                   "variant_regions": args.bed_file,
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
        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)

    workdir = utils.safe_makedir(os.path.join(args.base_dir, args.sample, "work"))
    return workdir, {"config_file": args.bcbio_config_file,
                     "run_info_yaml": out_config_file}
