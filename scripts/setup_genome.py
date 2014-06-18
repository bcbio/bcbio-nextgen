"""
Script to set up a custom genome for bcbio-nextgen
"""

import argparse
from argparse import ArgumentParser
import os
from cloudbio.biodata import genomes
from cloudbio import fabutils
from bcbio.utils import safe_makedir
from bcbio.distributed.transaction import file_transaction
from fabric.api import *
import subprocess
import sys
import shutils

SEQ_DIR = "seq"
RNASEQ_DIR = "rnaseq"

def _index_w_command(dir_name, command, ref_file, ext=None):
    index_name = os.path.splitext(os.path.basename(ref_file))[0]
    if ext is not None: index_name += ext
    build_path = os.path.join(os.path.dirname(ref_file), os.pardir)
    out_dir = os.path.join(build_path, dir_name)
    index_path = os.path.join(out_dir, index_name)
    if not env.safe_exists(out_dir):
        env.safe_run("mkdir %s" % out_dir)
    subprocess.check_call(command.format(ref_file=ref_file,
                                         index_name=index_path), shell=True)
    return index_path

genomes._index_w_command = _index_w_command

def setup_base_directories(genome_dir, name, build, gtf=None):
    name_dir = os.path.join(genome_dir, name)
    safe_makedir(name_dir)
    build_dir = os.path.join(name_dir, build)
    safe_makedir(build_dir)
    seq_dir = os.path.join(build_dir, SEQ_DIR)
    safe_makedir(seq_dir)
    if gtf:
        gtf_dir = os.path.join(build_dir, RNASEQ_DIR)
        safe_makedir(gtf_dir)
    return build_dir

def install_fasta_file(build_dir, fasta, build):
    out_file = os.path.join(build_dir, SEQ_DIR, build + ".fa")
    if not os.path.exists(out_file):
        shutils.copyfile(fasta, out_file)
    return out_file

def install_gtf_file(build_dir, gtf, build):
    out_file = os.path.join(build_dir, RNASEQ_DIR, "ref-transcripts.gtf")
    if not os.path.exists(out_file):
        shutils.copyfile(gtf, out_file)
    return out_file

def make_indices(ref_file, indices):
    return False

if __name__ == "__main__":
    description = ("Set up a custom genome for bcbio-nextgen. This will "
                   "stick the genome under name/build in the genomes "
                   "directory in your bcbio-nextgen installation.")

    parser = ArgumentParser(description=("Set up a custom genome for "
                                         "bcbio-nextgen."))
    parser.add_argument("-f", "--fasta", required=True,
                        help="FASTA file of the genome.")
    parser.add_argument("-g", "--gtf", default=None,
                        help="GTF file of the transcriptome")
    parser.add_argument("-n", "--name", required=True,
                        help="Name of genome.")
    parser.add_argument("-b", "--build", required=True,
                        help="Build of genome.")
    parser.add_argument("-d", "--genome-dir", required=True,
                        help="Path to bcbio-nextgen genomes directory.")
    parser.add_argument("-i", "--indexes", choices=genomes.INDEX_FNS.keys(),
                        nargs="*", required=True,
                        help="List of indexes to make")
    parser.add_argument("-p", "--picard-dir", default=None,
                        help="Path to Picard installation")
    parser.add_argument("--prepare-tx", default=None,
                        help="Path to prepare_tx_gff.py (in utils/cloudbiolinux)")

    args = parser.parse_args()
    fabutils.configure_runsudo(env)
    args.genome_dir = os.path.abspath(args.genome_dir)
    args.fasta = os.path.abspath(args.fasta)
    args.gtf = os.path.abspath(args.gtf) if args.gtf else None
    env.system_install = args.genome_dir

    print "Creating directories using %s as the base." % (args.genome_dir)
    build_dir = setup_base_directories(args.genome_dir, args.name,
                                       args.build, args.gtf)
    os.chdir(build_dir)
    print "Genomes will be installed into %s." % (build_dir)

    fasta_file = install_fasta_file(build_dir, args.fasta, args.build)
    print "Installed genome as %s." % (fasta_file)
    if args.gtf:
        gtf_file = install_gtf_file(build_dir, args.gtf, args.build)
        print "Installed GTF as %s." % (gtf_file)

    for index in args.indexes:
        print "Creating the %s index." % (index)
        index_fn = genomes.get_index_fn(index)
        if not index_fn:
            print "Do not know how to make the index %s, skipping." % (index)
        index_fn(fasta_file)

    if args.gtf:
        if not args.picard_dir or not args.prepare_tx:
            print ("picard-dir and prepare-tx must be set when preparing a "
                   "transcriptome.")
            sys.exit(1)

        "Preparing transcriptome."
        os.chdir(os.path.join(build_dir, os.pardir))
        cmd = ("python {args.prepare_tx} --gtf {gtf_file} {args.picard_dir} "
               "{args.build}")
        subprocess.check_call(cmd.format(**locals()), shell=True)
