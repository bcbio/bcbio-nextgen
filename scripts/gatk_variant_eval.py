#!/usr/bin/env python
"""Evaluate variant calls using GATK tools.

Usage:
    gatk_variant_eval.py <Picard location> <VCF file or directory> <reference seq file>
                         <dbSNP file> [<interval file of targets>]

If a directory is passed, this will print out a text table of SNP counts and
Transition/Transversion ratios. If a single file is passed, a JSON string will
be printed that could be picked up and parsed by calling programs in an
automated pipeline.

http://www.broadinstitute.org/gsa/wiki/index.php/VariantEval
"""
import os
import sys
import subprocess
import itertools
import glob
import json

import yaml

from bcbio.picard import PicardRunner

def main(picard_dir, vcf_info, ref_file, dbsnp, intervals=None):
    picard = PicardRunner(picard_dir)
    if os.path.isdir(vcf_info):
        vcf_files = sorted(glob.glob(os.path.join(vcf_info, "*-filter.vcf")))
    else:
        vcf_files = [vcf_info]
    for vcf_in in vcf_files:
        eval_file = variant_eval(vcf_in, ref_file, dbsnp, intervals, picard)
        stats = extract_eval_stats(eval_file)
        print_stats(vcf_in, stats['called'], len(vcf_files) == 1)

def print_stats(in_file, stats, print_json=False):
    in_info = "_".join(os.path.basename(in_file).split("_")[:2])
    total = sum(itertools.chain.from_iterable(s.itervalues() for s in stats.itervalues()))
    dbsnp = sum(stats['known'].itervalues()) / float(total) * 100.0
    tv_dbsnp = stats['known']['tv']
    ti_dbsnp = stats['known']['ti']
    tv_novel = stats['novel']['tv']
    ti_novel = stats['novel']['ti']
    titv_all = float(ti_novel + ti_dbsnp) / float(tv_novel + tv_dbsnp)
    titv_dbsnp = float(ti_dbsnp) / float(tv_dbsnp)
    titv_novel = float(ti_novel) / float(tv_novel)

    info = dict(total=total, dbsnp_pct = dbsnp, titv_all=titv_all,
            titv_dbsnp=titv_dbsnp, titv_novel=titv_novel)
    if print_json:
        print json.dumps(info)
    else:
        print "%s % 7s   %.1f    %.2f   %.2f   %.2f" % (in_info, total, dbsnp,
                titv_all, titv_dbsnp, titv_novel)

def extract_eval_stats(eval_file):
    """Parse statistics of interest from GATK output file.
    """
    stats = dict()
    for snp_type in ['called', 'filtered']:
        stats[snp_type]  = dict()
        for dbsnp_type in ['known', 'novel']:
            stats[snp_type][dbsnp_type] = dict(ti=0, tv=0)
    for line in _eval_analysis_type(eval_file, "Ti/Tv Variant Evaluator"):
        if line[:2] == ['eval', 'dbsnp']:
            snp_type = line[3]
            dbsnp_type = line[4]
            try:
                cur = stats[snp_type][dbsnp_type]
            except KeyError:
                cur = None
            if cur:
                stats[snp_type][dbsnp_type]["ti"] = int(line[5])
                stats[snp_type][dbsnp_type]["tv"] = int(line[6])
    return stats

def _eval_analysis_type(in_file, analysis_name):
    """Retrieve data lines associated with a particular analysis.
    """
    with open(in_file) as in_handle:
        # read until we reach the analysis
        for line in in_handle:
            if (line.startswith("Analysis Name:") and
                line.find(analysis_name) > 0):
                break
        # read off header lines
        for header in range(4):
            in_handle.next()
        # read the table until a blank line
        for line in in_handle:
            if not line.strip():
                break
            parts = line.rstrip("\n\r").split()
            yield parts

def variant_eval(vcf_in, ref_file, dbsnp, target_intervals, picard):
    """Evaluate variants in comparison with dbSNP reference.
    """
    out_file = "%s.eval" % os.path.splitext(vcf_in)[0]
    params = ["-T", "VariantEval",
              "-R", ref_file,
              "-B", "eval,VCF,%s" % vcf_in,
              "-B", "dbsnp,VCF,%s" % dbsnp,
              "-o", out_file,
              "-l", "INFO"
              ]
    if target_intervals:
        params.extend(["-L", target_intervals])
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        picard.run_gatk(params)
    return out_file

if __name__ == "__main__":
    main(*sys.argv[1:])
