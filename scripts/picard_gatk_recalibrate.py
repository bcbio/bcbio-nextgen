#!/usr/bin/env python
"""Provide GATK recalibration similarly to how it is done at Broad.

http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration
http://www.broadinstitute.org/gsa/wiki/index.php/Built-in_command-line_arguments
http://www.broadinstitute.org/gsa/wiki/index.php/The_DBSNP_rod

Usage:
    picard_maq_recalibrate.py <config YAML> <reference file> <align BAM file>
                              <snp file>

Process description, from Broad:
- Mark Duplicates

- Count Covariates:

java -Xmx4g -jar GATK-Picard.jar 
-T CountCovariates 
-cov ReadGroupCovariate 
-cov QualityScoreCovariate 
-cov CycleCovariate 
-cov DinucCovariate
-cov TileCovariate
-recalFile /path/to/analysis_dir/FLOWCELL.LANE.recal_data.csv 
-I /path/to/analysis_dir/FLOWCELL.LANE.aligned.duplicates_marked.bam 
-R /path/to/reference.fasta 
-l INFO 
--use_original_quals 
-U
-B dbsnp,PicardDbSNP,/path/to/reference.dbsnp

- GATK recalibration:

java -jar GATK-Picard.jar 
-T TableRecalibration 
-recalFile /path/to/analysis_dir/FLOWCELL.LANE.recal_data.csv 
-R /path/to/reference.fasta 
-I /path/to/analysis_dir/FLOWCELL.LANE.aligned.duplicates_marked.bam 
-outputBam /path/to/analysis_dir/FLOWCELL.LANE.gatk.bam 
-l INFO
--use_original_quals 
-U
"""
import os
import sys
import glob
import shutil
import subprocess

import yaml

from bcbio.broad import BroadRunner
from bcbio.utils import curdir_tmpdir

def main(config_file, ref_file, align_bam, snp_file=None):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    picard = BroadRunner(config["program"]["picard"], config["program"].get("gatk", ""))
    platform = config["algorithm"]["platform"]
    ref_dict = index_ref_file(picard, ref_file)
    #snp_dict = (index_snp_file(picard, ref_dict, snp_file) if snp_file else
    #        None)
    align_sort_bam = picard_sort(picard, align_bam)
    dup_align_bam = mark_duplicates(picard, align_sort_bam)
    recal_file = count_covariates(picard, dup_align_bam, ref_file, platform,
            snp_file)
    recal_bam = gatk_recalibrate(picard, dup_align_bam, ref_file, recal_file,
                                 platform)

def gatk_recalibrate(picard, dup_align_bam, ref_file, recal_file, platform):
    """Step 2 of GATK recalibration -- use covariates to re-write output file.
    """
    out_file = "%s-gatkrecal.bam" % os.path.splitext(dup_align_bam)[0]
    params = ["-T", "TableRecalibration",
              "-recalFile", recal_file,
              "-R", ref_file,
              "-I", dup_align_bam,
              "--out", out_file,
              "-baq",  "RECALCULATE",
              "-l", "INFO",
              "-U",
              "-OQ",
              "--default_platform", platform,
              ]
    if not os.path.exists(out_file):
        if _recal_available(recal_file):
            with curdir_tmpdir() as tmp_dir:
                picard.run_gatk(params, tmp_dir)
        else:
            shutil.copy(dup_align_bam, out_file)
    return out_file

def _recal_available(recal_file):
    """Determine if it's possible to do a recalibration; do we have data?
    """
    if os.path.exists(recal_file):
        with open(recal_file) as in_handle:
            while 1:
                line = in_handle.next()
                if not line.startswith("#"):
                    break
            test_line = in_handle.next()
            if test_line and not test_line.startswith("EOF"):
                return True
    return False

def count_covariates(picard, dup_align_bam, ref_file, platform,
        snp_file):
    """Step 1 of GATK recalibration process -- counting covariates.
    """
    out_file = "%s.recal" % os.path.splitext(dup_align_bam)[0]
    params = ["-T", "CountCovariates",
              "-cov", "ReadGroupCovariate",
              "-cov", "QualityScoreCovariate",
              "-cov", "CycleCovariate",
              "-cov", "DinucCovariate",
              "-cov", "TileCovariate",
              "-recalFile", out_file,
              "-I", dup_align_bam,
              "-R", ref_file,
              "-l", "INFO",
              "-U",
              "-OQ",
              "--default_platform", platform,
              ]
    if snp_file:
        params += ["-B:dbsnp,VCF", snp_file]
    if not os.path.exists(out_file):
        with curdir_tmpdir() as tmp_dir:
            picard.run_gatk(params, tmp_dir)
    return out_file

def mark_duplicates(picard, align_bam):
    base, ext = os.path.splitext(align_bam)
    base = base.replace(".", "-")
    dup_bam = "%s-dup%s" % (base, ext)
    dup_metrics = "%s-dup.dup_metrics" % base
    if not os.path.exists(dup_bam):
        with curdir_tmpdir() as tmp_dir:
            opts = [("INPUT", align_bam),
                    ("OUTPUT", dup_bam),
                    ("TMP_DIR", tmp_dir),
                    ("METRICS_FILE", dup_metrics)]
        picard.run("MarkDuplicates", opts)
    return dup_bam

def picard_sort(picard, align_bam):
    base, ext = os.path.splitext(align_bam)
    out_file = "%s-sort%s" % (base, ext)
    if not os.path.exists(out_file):
        with curdir_tmpdir() as tmp_dir:
            opts = [("INPUT", align_bam),
                    ("OUTPUT", out_file),
                    ("TMP_DIR", tmp_dir),
                    ("SORT_ORDER", "coordinate")]
            picard.run("SortSam", opts)
    return out_file

def NOTUSED_index_snp_file(picard, ref_dict, snp_file):
    """Provide a DbSNP index file to pass to downstream analyses.

    XXX Drop this as soon as we are happy with VCF formatted dbSNP.
    """
    snp_index = "%s.dbsnp" % os.path.splitext(snp_file)[0]
    if not os.path.exists(snp_index):
        opts = [("SNP_FILE", snp_file),
                ("SEQUENCE_DICTIONARY", ref_dict),
                ("OUTPUT", snp_index)]
        picard.run("GenerateDbSnpFile", opts)
    return snp_index

def index_ref_file(picard, ref_file):
    """Provide a Picard style dict index file for a reference genome.
    """
    dict_file = "%s.dict" % os.path.splitext(ref_file)[0]
    if not os.path.exists(dict_file):
        opts = [("REFERENCE", ref_file),
                ("OUTPUT", dict_file)]
        picard.run("CreateSequenceDictionary", opts)
    return dict_file

if __name__ == "__main__":
    main(*sys.argv[1:])
