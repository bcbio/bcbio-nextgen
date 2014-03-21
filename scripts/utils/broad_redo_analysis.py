#!/usr/bin/env python
"""Redo post-processing of Broad alignments with updated pipeline.

Usage:
    broad_redo_analysis.py <YAML config file> <flow cell dir>
"""
import os
import sys
import json
import contextlib
import subprocess
import glob
import copy
import csv
from optparse import OptionParser
from multiprocessing import Pool
import xml.etree.ElementTree as ET

import yaml

from bcbio.illumina import flowcell
from bcbio.galaxy.api import GalaxyApiAccess
from bcbio.broad.metrics import PicardMetricsParser
from bcbio import utils
from bcbio.pipeline.config_utils import load_config

def main(config_file, fc_dir):
    work_dir = os.getcwd()
    config = load_config(config_file)
    galaxy_api = GalaxyApiAccess(config['galaxy_url'], config['galaxy_api_key'])
    fc_name, fc_date = flowcell.parse_dirname(fc_dir)
    run_info = galaxy_api.run_details(fc_name)
    fastq_dir = flowcell.get_fastq_dir(fc_dir)
    if config["algorithm"]["num_cores"] > 1:
        pool = Pool(config["algorithm"]["num_cores"])
        try:
            pool.map(_process_wrapper,
                    ((i, fastq_dir, fc_name, fc_date, config, config_file)
                        for i in run_info["details"]))
        except:
            pool.terminate()
            raise
    else:
        map(_process_wrapper,
            ((i, fastq_dir, fc_name, fc_date, config, config_file)
                for i in run_info["details"]))

def process_lane(info, fastq_dir, fc_name, fc_date, config, config_file):
    config = _update_config_w_custom(config, info)
    sample_name = info.get("description", "")
    if config["algorithm"]["include_short_name"]:
        sample_name = "%s: %s" % (info.get("name", ""), sample_name)
    genome_build = "%s%s" % (info["genome_build"],
                             config["algorithm"].get("ref_ext", ""))
    if info.get("analysis", "") == "Broad SNP":
        print "Processing", info["lane"], genome_build, \
                sample_name, info.get("researcher", ""), \
                info.get("analysis", "")
        lane_name = "%s_%s_%s" % (info['lane'], fc_date, fc_name)
        aligner_to_use = config["algorithm"]["aligner"]
        align_ref, sam_ref = get_genome_ref(genome_build,
                aligner_to_use, os.path.dirname(config["galaxy_config"]))
        util_script_dir = os.path.dirname(__file__)
        base_bam = "%s.bam" % lane_name
        resort_bam = resort_karotype_and_rename(base_bam, sam_ref, util_script_dir)
        base_bam = resort_bam
        if config["algorithm"]["recalibrate"]:
            print info['lane'], "Recalibrating with GATK"
            dbsnp_file = get_dbsnp_file(config, sam_ref)
            gatk_bam = recalibrate_quality(base_bam, sam_ref,
                    dbsnp_file, config["program"]["picard"])
            if config["algorithm"]["snpcall"]:
                print info['lane'], "Providing SNP genotyping with GATK"
                run_genotyper(gatk_bam, sam_ref, dbsnp_file, config_file)

def resort_karotype_and_rename(in_bam, ref_file, script_dir):
    assert os.path.exists(in_bam)
    out_file = "%s-ksort%s" % os.path.splitext(in_bam)
    ref_dict = "%s.dict" % os.path.splitext(ref_file)[0]
    assert os.path.exists(ref_dict)
    if not os.path.exists(out_file):
        resort_script = os.path.join(script_dir, "resort_bam_karyotype.py")
        rename_script = os.path.join(script_dir, "rename_samples.py")
        print "Resorting to karyotype", in_bam
        cl = ["python2.6", resort_script, ref_dict, in_bam]
        subprocess.check_call(cl)
        print "Renaming samples", out_file
        cl = ["python2.6", rename_script, out_file]
        subprocess.check_call(cl)
    return out_file

def _process_wrapper(args):
    try:
        return process_lane(*args)
    except KeyboardInterrupt:
        raise Exception

def recalibrate_quality(bam_file, sam_ref, dbsnp_file, picard_dir):
    """Recalibrate alignments with GATK and provide pdf summary.
    """
    cl = ["picard_gatk_recalibrate.py", picard_dir, sam_ref, bam_file]
    if dbsnp_file:
        cl.append(dbsnp_file)
    subprocess.check_call(cl)
    out_file = glob.glob("%s*gatkrecal.bam" % os.path.splitext(bam_file)[0])[0]
    return out_file

def run_genotyper(bam_file, ref_file, dbsnp_file, config_file):
    """Perform SNP genotyping and analysis using GATK.
    """
    cl = ["gatk_genotyper.py", config_file, ref_file, bam_file]
    if dbsnp_file:
        cl.append(dbsnp_file)
    subprocess.check_call(cl)

def get_dbsnp_file(config, sam_ref):
    snp_file = config["algorithm"].get("dbsnp", None)
    if snp_file:
        base_dir = os.path.dirname(os.path.dirname(sam_ref))
        snp_file = os.path.join(base_dir, snp_file)
    return snp_file

def get_fastq_files(directory, lane, fc_name):
    """Retrieve fastq files for the given lane, ready to process.
    """
    files = glob.glob(os.path.join(directory, "%s_*%s*txt*" % (lane, fc_name)))
    files.sort()
    if len(files) > 2 or len(files) == 0:
        raise ValueError("Did not find correct files for %s %s %s %s" %
                (directory, lane, fc_name, files))
    ready_files = []
    for fname in files:
        if fname.endswith(".gz"):
            cl = ["gunzip", fname]
            subprocess.check_call(cl)
            ready_files.append(os.path.splitext(fname)[0])
        else:
            ready_files.append(fname)
    return ready_files[0], (ready_files[1] if len(ready_files) > 1 else None)

def _remap_to_maq(ref_file):
    base_dir = os.path.dirname(os.path.dirname(ref_file))
    name = os.path.basename(ref_file)
    for ext in ["fa", "fasta"]:
        test_file = os.path.join(base_dir, "maq", "%s.%s" % (name, ext))
        if os.path.exists(test_file):
            return test_file
    raise ValueError("Did not find maq file %s" % ref_file)

def get_genome_ref(genome_build, aligner, galaxy_base):
    """Retrieve the reference genome file location from galaxy configuration.
    """
    ref_files = dict(
            bowtie = "bowtie_indices.loc",
            bwa = "bwa_index.loc",
            samtools = "sam_fa_indices.loc",
            maq = "bowtie_indices.loc")
    remap_fns = dict(
            maq = _remap_to_maq
            )
    out_info = []
    for ref_get in [aligner, "samtools"]:
        ref_file = os.path.join(galaxy_base, "tool-data", ref_files[ref_get])
        with open(ref_file) as in_handle:
            for line in in_handle:
                if not line.startswith("#"):
                    parts = line.strip().split()
                    if parts[0] == "index":
                        parts = parts[1:]
                    if parts[0] == genome_build:
                        out_info.append(parts[-1])
                        break
        try:
            out_info[-1] = remap_fns[ref_get](out_info[-1])
        except KeyError:
            pass
        except IndexError:
            raise IndexError("Genome %s not found in %s" % (genome_build,
                ref_file))

    if len(out_info) != 2:
        raise ValueError("Did not find genome reference for %s %s" %
                (genome_build, aligner))
    else:
        return tuple(out_info)

# Utility functions

def _update_config_w_custom(config, lane_info):
    """Update the configuration for this lane if a custom analysis is specified.
    """
    config = copy.deepcopy(config)
    custom = config["custom_algorithms"].get(lane_info.get("analysis", None),
            None)
    if custom:
        for key, val in custom.iteritems():
            config["algorithm"][key] = val
    return config

if __name__ == "__main__":
    parser = OptionParser()
    (options, args) = parser.parse_args()
    kwargs = dict()
    main(*args, **kwargs)
