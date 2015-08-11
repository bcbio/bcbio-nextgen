#!/usr/bin/env python -Es
"""
Script to set up a custom genome for bcbio-nextgen
"""

import argparse
from argparse import ArgumentParser
import os
import toolz as tz
from bcbio.utils import safe_makedir, file_exists, chdir
from bcbio.pipeline import config_utils
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.install import (REMOTES, get_cloudbiolinux, SUPPORTED_GENOMES, SUPPORTED_INDEXES,
                           _get_data_dir)
from bcbio.galaxy import loc
from fabric.api import *
import subprocess
import sys
import shutil
import yaml
import gffutils
from gffutils.iterators import DataIterator
import tempfile

SEQ_DIR = "seq"
RNASEQ_DIR = "rnaseq"
SRNASEQ_DIR = "srnaseq"

ERCC_BUCKET = "bcbio-data.s3.amazonaws.com/"

def gff3_to_gtf(gff3_file):

    dialect = {'field separator': '; ',
               'fmt': 'gtf',
               'keyval separator': ' ',
               'leading semicolon': False,
               'multival separator': ',',
               'quoted GFF2 values': True,
               'order': ['gene_id', 'transcript_id'],
               'repeated keys': False,
               'trailing semicolon': True}

    out_file = os.path.splitext(gff3_file)[0] + ".gtf"
    if file_exists(out_file):
        return out_file

    print "Converting %s to %s." %(gff3_file, out_file)

    db = gffutils.create_db(gff3_file, ":memory:")
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            for feature in DataIterator(db.features_of_type("exon"), dialect=dialect):
                transcript_id = feature["Parent"][0]
                gene_id = db[transcript_id]["Parent"][0]
                attr = {"transcript_id": transcript_id, "gene_id": gene_id}
                attributes = gffutils.attributes.Attributes(attr)
                feature.attributes = attributes
                print >> out_handle, feature
    return out_file


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
        shutil.copyfile(fasta, out_file)
    return out_file

def install_gtf_file(build_dir, gtf, build):
    out_file = os.path.join(build_dir, RNASEQ_DIR, "ref-transcripts.gtf")
    if not os.path.exists(out_file):
        shutil.copyfile(gtf, out_file)
    return out_file

def install_srna(species, gtf):
    out_file = os.path.join(SRNASEQ_DIR, "srna-transcripts.gtf")
    safe_makedir(SRNASEQ_DIR)
    if not os.path.exists(out_file):
        shutil.copyfile(gtf, out_file)
    try:
        from seqcluster import install
    except ImportError:
        raise ImportError("install seqcluster first, please.")
    with chdir(SRNASEQ_DIR):
        hairpin, miRNA = install._install_mirbase()
        cmd = ("grep -A 2 {species} {hairpin} | grep -v '\-\-$' | tr U T  > hairpin.fa")
        do.run(cmd.format(**locals()), "set precursor.")
        cmd = ("grep -A 1 {species} {miRNA} > miRNA.str")
        do.run(cmd.format(**locals()), "set miRNA.")
        shutil.rmtree("mirbase")
    return out_file

def append_ercc(gtf_file, fasta_file):
    ercc_fa = ERCC_BUCKET + "ERCC92.fasta.gz"
    tmp_fa = tempfile.NamedTemporaryFile(delete=False, suffix=".gz").name
    append_fa_cmd = "wget {ercc_fa} -O {tmp_fa}; gzip -cd {tmp_fa} >> {fasta_file}"
    print append_fa_cmd.format(**locals())
    subprocess.check_call(append_fa_cmd.format(**locals()), shell=True)
    ercc_gtf = ERCC_BUCKET + "ERCC92.gtf.gz"
    tmp_gtf = tempfile.NamedTemporaryFile(delete=False, suffix=".gz").name
    append_gtf_cmd = "wget {ercc_gtf} -O {tmp_gtf}; gzip -cd {tmp_gtf} >> {gtf_file}"
    print append_gtf_cmd.format(**locals())
    subprocess.check_call(append_gtf_cmd.format(**locals()), shell=True)


if __name__ == "__main__":
    description = ("Set up a custom genome for bcbio-nextgen. This will "
                   "place the genome under name/build in the genomes "
                   "directory in your bcbio-nextgen installation.")

    parser = ArgumentParser(description=description)
    parser.add_argument("-f", "--fasta", required=True,
                        help="FASTA file of the genome.")
    parser.add_argument("--gff3", default=False, action='store_true',
                        help="File is a GFF3 file.")
    parser.add_argument("-g", "--gtf", default=None,
                        help="GTF file of the transcriptome")
    parser.add_argument("-n", "--name", required=True,
                        help="Name of organism, for example Hsapiens.")
    parser.add_argument("-b", "--build", required=True,
                        help="Build of genome, for example hg19.")
    parser.add_argument("-i", "--indexes", choices=SUPPORTED_INDEXES, nargs="*",
                        default=["seq"], help="Space separated list of indexes to make")
    parser.add_argument("--ercc", action='store_true', default=False,
                        help="Add ERCC spike-ins.")
    parser.add_argument("--mirbase", help="species in mirbase for smallRNAseq data.")
    parser.add_argument("--srna_gtf", help="gtf to use for smallRNAseq data.")

    args = parser.parse_args()
    if not all([args.mirbase, args.srna_gtf]) and any([args.mirbase, args.srna_gtf]):
        raise ValueError("--mirbase and --srna_gtf both need a value.")

    env.hosts = ["localhost"]
    cbl = get_cloudbiolinux(REMOTES)
    sys.path.insert(0, cbl["dir"])
    genomemod = __import__("cloudbio.biodata", fromlist=["genomes"])
    # monkey patch cloudbiolinux to use this indexing command instead
    genomes = getattr(genomemod, 'genomes')
    genomes._index_w_command = _index_w_command
    fabmod = __import__("cloudbio", fromlist=["fabutils"])
    fabutils = getattr(fabmod, 'fabutils')
    fabutils.configure_runsudo(env)

    system_config = os.path.join(_get_data_dir(), "galaxy", "bcbio_system.yaml")
    with open(system_config) as in_handle:
        config = yaml.load(in_handle)
    env.picard_home = config_utils.get_program("picard", config, ptype="dir")

    genome_dir = os.path.abspath(os.path.join(_get_data_dir(), "genomes"))
    args.fasta = os.path.abspath(args.fasta)
    args.gtf = os.path.abspath(args.gtf) if args.gtf else None
    if args.gff3:
        args.gtf = gff3_to_gtf(args.gtf)

    # always make a sequence dictionary
    if "seq" not in args.indexes:
        args.indexes.append("seq")

    env.system_install = genome_dir
    prepare_tx = os.path.join(cbl["dir"], "utils", "prepare_tx_gff.py")

    print "Creating directories using %s as the base." % (genome_dir)
    build_dir = setup_base_directories(genome_dir, args.name, args.build, args.gtf)
    os.chdir(build_dir)
    print "Genomes will be installed into %s." % (build_dir)

    fasta_file = install_fasta_file(build_dir, args.fasta, args.build)
    print "Installed genome as %s." % (fasta_file)
    if args.gtf:
        if "bowtie2" not in args.indexes:
            args.indexes.append("bowtie2")
        gtf_file = install_gtf_file(build_dir, args.gtf, args.build)
        print "Installed GTF as %s." % (gtf_file)

    if args.ercc:
        print "Appending ERCC sequences to %s and %s." % (gtf_file, fasta_file)
        append_ercc(gtf_file, fasta_file)

    indexed = {}
    for index in args.indexes:
        print "Creating the %s index." % (index)
        index_fn = genomes.get_index_fn(index)
        if not index_fn:
            print "Do not know how to make the index %s, skipping." % (index)
            continue
        indexed[index] = index_fn(fasta_file)
    indexed["samtools"] = fasta_file

    if args.gtf:
        "Preparing transcriptome."
        with chdir(os.path.join(build_dir, os.pardir)):
            cmd = ("{sys.executable} {prepare_tx} --gtf {gtf_file} {args.build}")
            subprocess.check_call(cmd.format(**locals()), shell=True)
    if args.mirbase:
        "Preparing smallRNA data."
        with chdir(os.path.join(build_dir)):
            install_srna(args.mirbase, args.srna_gtf)

    base_dir = os.path.normpath(os.path.dirname(fasta_file))
    resource_file = os.path.join(base_dir, "%s-resources.yaml" % args.build)

    print "Dumping genome resources to %s." % resource_file
    resource_dict = {"version": 1}
    if args.gtf:
        transcripts = ["rnaseq", "transcripts"]
        mask = ["rnaseq", "transcripts_mask"]
        index = ["rnaseq", "transcriptome_index", "tophat"]
        dexseq = ["rnaseq", "dexseq"]
        refflat = ["rnaseq", "refflat"]
        rRNA_fa = ["rnaseq", "rRNA_fa"]
        resource_dict = tz.update_in(resource_dict, transcripts,
                                     lambda x: "../rnaseq/ref-transcripts.gtf")
        resource_dict = tz.update_in(resource_dict, mask,
                                     lambda x: "../rnaseq/ref-transcripts-mask.gtf")
        resource_dict = tz.update_in(resource_dict, index,
                                     lambda x: "../rnaseq/tophat/%s_transcriptome.ver" % args.build)
        resource_dict = tz.update_in(resource_dict, refflat,
                                     lambda x: "../rnaseq/ref-transcripts.refFlat")
        resource_dict = tz.update_in(resource_dict, dexseq,
                                     lambda x: "../rnaseq/ref-transcripts.dexseq.gff3")
        resource_dict = tz.update_in(resource_dict, rRNA_fa,
                                     lambda x: "../rnaseq/rRNA.fa")
    if args.mirbase:
        srna_gtf = ["srnaseq", "srna-transcripts"]
        srna_mirbase = ["srnaseq", "mirbase"]
        resource_dict = tz.update_in(resource_dict, srna_gtf,
                                     lambda x: "../srnaseq/srna-transcripts.gtf")
        resource_dict = tz.update_in(resource_dict, srna_mirbase,
                                     lambda x: "../srnaseq/hairpin.fa")
    # write out resource dictionarry
    with file_transaction(resource_file) as tx_resource_file:
        with open(tx_resource_file, "w") as out_handle:
            out_handle.write(yaml.dump(resource_dict, default_flow_style=False))

    print "Updating Galaxy .loc files."
    galaxy_base = os.path.join(_get_data_dir(), "galaxy")
    for index, index_file in indexed.items():
        loc.update_loc_file(galaxy_base, index, args.build, index_file)
