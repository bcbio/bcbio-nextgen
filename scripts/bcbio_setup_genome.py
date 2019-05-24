#!/usr/bin/env python -Es
"""
Script to set up a custom genome for bcbio-nextgen
"""
from __future__ import print_function

from argparse import ArgumentParser
import collections
import gzip
import os
from Bio import SeqIO
import toolz as tz

from bcbio.utils import safe_makedir, file_exists, chdir, is_gzipped
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do

from bcbio.install import (REMOTES, get_cloudbiolinux, SUPPORTED_INDEXES,
                           _get_data_dir)
from bcbio.pipeline.run_info import ALLOWED_CONTIG_NAME_CHARS
from bcbio.galaxy import loc
from bcbio.log import logger

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

def extract_if_gzipped(filename):
    stem, ext = os.path.splitext(filename)
    if ext == ".gz":
        subprocess.check_call("gzip -cd %s > %s" % (filename, stem), shell=True)
        return stem
    else:
        return filename

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

    logger.info("Converting %s to %s." % (gff3_file, out_file))

    if _is_from_ncbi(gff3_file):
        logger.info("NCBI format detected by the presence of the %s key."
                    % _is_from_ncbi(gff3_file))
        _output_ncbi_gff3(gff3_file, out_file, dialect)
    else:
        _output_gff3(gff3_file, out_file, dialect)
    return out_file

def _output_gff3(gff3_file, out_file, dialect):
    db = gffutils.create_db(gff3_file, ":memory:")
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            for feature in DataIterator(db.features_of_type("exon"), dialect=dialect):
                transcript_id = feature["Parent"][0]
                gene_id = db[transcript_id]["Parent"][0]
                attr = {"transcript_id": transcript_id, "gene_id": gene_id}
                attributes = gffutils.attributes.Attributes(attr)
                feature.attributes = attributes
                print(feature, file=out_handle, end="")

def _output_ncbi_gff3(gff3_file, out_file, dialect):
    gene_key = "gene"
    id_spec = {"gene": gene_key}
    db = gffutils.create_db(gff3_file, ":memory:", id_spec=id_spec)
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            for feature in DataIterator(db.features_of_type("exon"), dialect=dialect):
                # Gnomon features are often missing a transcript id
                # some malformed features are also missing the gene key
                try:
                    transcript_id = feature["transcript_id"]
                except KeyError:
                    try:
                        transcript_id = feature[gene_key]
                    except KeyError:
                        continue
                gene_id = feature[gene_key]
                try:
                    biotype = feature["gene_biotype"]
                except KeyError:
                    biotype = "unknown"
                attr = {"transcript_id": transcript_id, "gene_id": gene_id,
                        "gene_biotype": biotype}
                attributes = gffutils.attributes.Attributes(attr)
                feature.attributes = attributes
                print(feature, file=out_handle, end="")

def _is_from_ncbi(gff3_file):
    with open(gff3_file) as in_handle:
        for line in tz.take(10000, in_handle):
            if "Dbxref" in line:
                return "Dbxref"
            if "db_xref" in line:
                return "db_xref"
    return None

def _index_w_command(env, dir_name, command, ref_file, pre=None, post=None, ext=None):
    index_name = os.path.splitext(os.path.basename(ref_file))[0]
    if ext is not None: index_name += ext
    build_path = os.path.join(os.path.dirname(ref_file), os.pardir)
    out_dir = os.path.join(build_path, dir_name)
    index_path = os.path.join(out_dir, index_name)
    safe_makedir(out_dir)
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
    if not file_exists(out_file):
        recs = SeqIO.parse(fasta, "fasta")
        with open(out_file, "w") as out_handle:
            SeqIO.write((_clean_rec_name(rec) for rec in recs), out_handle, "fasta")
    return out_file

def _clean_rec_name(rec):
    """Clean illegal characters in input fasta file which cause problems downstream.
    """
    out_id = []
    for char in list(rec.id):
        if char in ALLOWED_CONTIG_NAME_CHARS:
            out_id.append(char)
        else:
            out_id.append("_")
    rec.id = "".join(out_id)
    rec.description = ""
    return rec

def install_gtf_file(build_dir, gtf, build):
    out_file = os.path.join(build_dir, RNASEQ_DIR, "ref-transcripts.gtf")
    if not file_exists(out_file):
        if is_gzipped(gtf):
            with gzip.open(gtf_file, 'rb') as in_handle:
                with open(out_file, 'wb') as out_handle:
                    shutil.copyfileobj(in_handle, out_handle)
        else:
            shutil.copyfile(gtf, out_file)
    return out_file

def install_srna(species, gtf):
    out_file = os.path.join(SRNASEQ_DIR, "srna-transcripts.gtf")
    safe_makedir(SRNASEQ_DIR)
    if gtf:
        if not file_exists(out_file):
            shutil.copyfile(gtf, out_file)
    try:
        from seqcluster import install
    except ImportError:
        raise ImportError("install seqcluster first, please.")
    with chdir(SRNASEQ_DIR):
        hairpin, miRNA = install._install_mirbase()
        cmd = ("cat %s |  awk '{if ($0~/>%s/){name=$0; print name} else if ($0~/^>/){name=0};if (name!=0 && $0!~/^>/){print $0;}}' | sed 's/U/T/g'  > hairpin.fa")
        do.run(cmd % (hairpin, species), "set precursor.")
        cmd = ("grep -A 1 {species} {miRNA} > miRNA.str")
        do.run(cmd.format(**locals()), "set miRNA.")
        shutil.rmtree("mirbase")
    return out_file

def append_ercc(gtf_file, fasta_file):
    ercc_fa = ERCC_BUCKET + "ERCC92.fasta.gz"
    tmp_fa = tempfile.NamedTemporaryFile(delete=False, suffix=".gz").name
    append_fa_cmd = "wget {ercc_fa} -O {tmp_fa}; gzip -cd {tmp_fa} >> {fasta_file}"
    print(append_fa_cmd.format(**locals()))
    subprocess.check_call(append_fa_cmd.format(**locals()), shell=True)
    ercc_gtf = ERCC_BUCKET + "ERCC92.gtf.gz"
    tmp_gtf = tempfile.NamedTemporaryFile(delete=False, suffix=".gz").name
    append_gtf_cmd = "wget {ercc_gtf} -O {tmp_gtf}; gzip -cd {tmp_gtf} >> {gtf_file}"
    print(append_gtf_cmd.format(**locals()))
    subprocess.check_call(append_gtf_cmd.format(**locals()), shell=True)

class MyParser(ArgumentParser):
    def error(self, message):
        self.print_help()
        galaxy_base = os.path.join(_get_data_dir(), "galaxy")
        print("\nCurrent genomes\n")
        print(open(loc.get_loc_file(galaxy_base, "samtools")).read())
        sys.exit(0)


if __name__ == "__main__":
    description = ("Set up a custom genome for bcbio-nextgen. This will "
                   "place the genome under name/build in the genomes "
                   "directory in your bcbio-nextgen installation.")

    parser = MyParser(description=description)

    parser.add_argument("-c", "--cores", default=1,
                        help="number of cores to use")
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
 #   if not all([args.mirbase, args.srna_gtf]) and any([args.mirbase, args.srna_gtf]):
 #       raise ValueError("--mirbase and --srna_gtf both need a value.")

    os.environ["PATH"] += os.pathsep + os.path.dirname(sys.executable)
    cbl = get_cloudbiolinux(REMOTES)
    sys.path.insert(0, cbl["dir"])
    genomemod = __import__("cloudbio.biodata", fromlist=["genomes"])
    # monkey patch cloudbiolinux to use this indexing command instead
    genomes = getattr(genomemod, 'genomes')
    genomes._index_w_command = _index_w_command

    genome_dir = os.path.abspath(os.path.join(_get_data_dir(), "genomes"))
    args.fasta = os.path.abspath(args.fasta)
    if not file_exists(args.fasta):
        print("%s does not exist, exiting." % args.fasta)
        sys.exit(1)

    args.gtf = os.path.abspath(args.gtf) if args.gtf else None
    if args.gtf and not file_exists(args.gtf):
        print("%s does not exist, exiting." % args.gtf)
        sys.exit(1)
    args.srna_gtf = os.path.abspath(args.srna_gtf) if args.srna_gtf else None

    gtf_file = args.gtf
    if args.gff3:
        gtf_file = extract_if_gzipped(gtf_file)
        gtf_file = gff3_to_gtf(gtf_file)

    # always make a sequence dictionary
    if "seq" not in args.indexes:
        args.indexes.append("seq")

    prepare_tx = os.path.join(cbl["dir"], "utils", "prepare_tx_gff.py")

    print("Creating directories using %s as the base." % (genome_dir))
    build_dir = setup_base_directories(genome_dir, args.name, args.build, args.gtf)
    os.chdir(build_dir)
    print("Genomes will be installed into %s." % (build_dir))

    fasta_file = extract_if_gzipped(args.fasta)
    fasta_file = install_fasta_file(build_dir, fasta_file, args.build)
    print("Installed genome as %s." % (fasta_file))
    if args.gtf:
        if "bowtie2" not in args.indexes:
            args.indexes.append("bowtie2")
        gtf_file = install_gtf_file(build_dir, gtf_file, args.build)
        print("Installed GTF as %s." % (gtf_file))

    if args.ercc:
        print("Appending ERCC sequences to %s and %s." % (gtf_file, fasta_file))
        append_ercc(gtf_file, fasta_file)

    indexed = {}
    Env = collections.namedtuple("Env", "system_install, cores")
    env = Env(genome_dir, args.cores)
    for index in args.indexes:
        print("Creating the %s index." % (index))
        index_fn = genomes.get_index_fn(index)
        if not index_fn:
            print("Do not know how to make the index %s, skipping." % (index))
            continue
        indexed[index] = index_fn(env, fasta_file)
    indexed["samtools"] = fasta_file

    if args.gtf:
        "Preparing transcriptome."
        with chdir(os.path.join(build_dir, os.pardir)):
            cmd = ("{sys.executable} {prepare_tx} --cores {args.cores} --genome-dir {genome_dir} "
                   "--gtf {gtf_file} {args.name} {args.build}")
            subprocess.check_call(cmd.format(**locals()), shell=True)
    if args.mirbase:
        "Preparing smallRNA data."
        with chdir(os.path.join(build_dir)):
            install_srna(args.mirbase, args.srna_gtf)

    base_dir = os.path.normpath(os.path.dirname(fasta_file))
    resource_file = os.path.join(base_dir, "%s-resources.yaml" % args.build)

    print("Dumping genome resources to %s." % resource_file)
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
        srna_gtf = ["srnaseq", "srna_transcripts"]
        srna_mirbase = ["srnaseq", "mirbase_hairpin"]
        resource_dict = tz.update_in(resource_dict, srna_gtf,
                                     lambda x: "../srnaseq/srna-transcripts.gtf")
        resource_dict = tz.update_in(resource_dict, srna_mirbase,
                                     lambda x: "../srnaseq/hairpin.fa")
    # write out resource dictionarry
    with file_transaction(resource_file) as tx_resource_file:
        with open(tx_resource_file, "w") as out_handle:
            out_handle.write(yaml.dump(resource_dict, default_flow_style=False))

    print("Updating Galaxy .loc files.")
    galaxy_base = os.path.join(_get_data_dir(), "galaxy")
    for index, index_file in indexed.items():
        if index_file:
            loc.update_loc_file(galaxy_base, index, args.build, index_file)

    print("Genome installation complete.")
