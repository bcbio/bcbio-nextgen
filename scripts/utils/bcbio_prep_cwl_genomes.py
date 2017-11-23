#!/usr/bin/env python
"""Clean and prepare a set of genomes for CWL usage and upload.

bcbio with CWL can read directly from a reference genome folder
without using Galaxy location files. This allows both local and
remote usage on object stores (Arvados, DNAnexus, SevenBridges, Synapse, S3).

This copies from an existing bcbio genome installation, cleaning
and packing directories to be ready for CWL usage and upload.

Usage:
      bcbio_prep_cwl_genomes.py <genome_dir>
"""
import glob
import os
import shutil
import subprocess
import sys
import tarfile

from bcbio import utils

def main(base_dir):
    for genome_dir in sorted(glob.glob(os.path.join(base_dir, "*", "*"))):
        if os.path.isdir(genome_dir):
            genome_name = os.path.basename(genome_dir)
            genome_out_dir = utils.safe_makedir(os.path.join(os.path.join(os.getcwd(), "genomes", genome_name)))
            copy_genome(genome_dir, genome_out_dir)

def copy_genome(orig_dir, out_dir):
    print(orig_dir, out_dir)
    to_copy = ["versions.csv", "bwa", "config", "coverage", "rnaseq", "rtg", "seq", "snpeff",
               "ucsc", "validation", "variation", "viral"]
    excludes = {"seq": ["*.fa.gz*", "*.old*", "perl"],
                "rnaseq": ["ericscript", "tophat", "kallisto"],
                "snpeff": ["transcripts"],
                "variation": ["genesplicer", "dbNSFP*"]}
    to_tar = ["bwa", "rtg", "snpeff"]
    for copy in to_copy:
        if os.path.isfile(os.path.join(orig_dir, copy)):
            shutil.copy(os.path.join(orig_dir, copy), out_dir)
        elif copy in to_tar and len(glob.glob(os.path.join(out_dir, "%s*-wf.tar.gz" % copy))) == 1:
            print("already prepped: %s" % glob.glob(os.path.join(out_dir, "%s*-wf.tar.gz" % copy)))
        else:
            cmd = ["rsync", "-avz"]
            for e in excludes.get(copy, []):
                cmd += ["--exclude", e]
            cmd += ["%s/%s/" % (orig_dir, copy), "%s/%s/" % (out_dir, copy)]
            print " ".join(cmd)
            subprocess.check_call(cmd)
            if copy in to_tar:
                with utils.chdir(out_dir):
                    out_file = copy
                    dir_files = os.listdir(copy)
                    if len(dir_files) == 1 and os.path.isdir(os.path.join(copy, dir_files[0])):
                        out_file += "--%s" % (dir_files[0])
                    out_file += "-wf.tar.gz"
                    print("tarball", out_file)
                    with tarfile.open(out_file, "w:gz") as tar:
                        tar.add(copy)
                    shutil.rmtree(copy)

if __name__ == "__main__":
    main(*sys.argv[1:])
