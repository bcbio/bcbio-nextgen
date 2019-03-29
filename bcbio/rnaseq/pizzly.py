"""
run the pizzly fusion caller for RNA-seq
https://github.com/pmelsted/pizzly
http://www.biorxiv.org/content/early/2017/07/20/166322
"""
from __future__ import print_function

import os

from bcbio.log import logger
from bcbio import utils
import bcbio.pipeline.datadict as dd
from bcbio.pipeline import config_utils
from bcbio.distributed.transaction import file_transaction
from bcbio.rnaseq import kallisto, sailfish, gtf
from bcbio.provenance import do
from bcbio.utils import file_exists, safe_makedir
from bcbio.bam import fasta

h5py = utils.LazyImport("h5py")
import numpy as np
import pandas as pd

def get_fragment_length(data):
    """
    lifted from
    https://github.com/pmelsted/pizzly/scripts/pizzly_get_fragment_length.py
    """
    h5 = kallisto.get_kallisto_h5(data)
    cutoff = 0.95
    with h5py.File(h5) as f:
        x = np.asarray(f['aux']['fld'], dtype='float64')
    y = np.cumsum(x)/np.sum(x)
    fraglen = np.argmax(y > cutoff)
    return(fraglen)

def run_pizzly(data):
    samplename = dd.get_sample_name(data)
    work_dir = dd.get_work_dir(data)
    pizzlydir = os.path.join(work_dir, "pizzly")
    gtf = dd.get_transcriptome_gtf(data)
    if not gtf:
        gtf = dd.get_gtf_file(data)
    if dd.get_transcriptome_fasta(data):
        gtf_fa = dd.get_transcriptome_fasta(data)
    else:
        gtf_fa = sailfish.create_combined_fasta(data)
    stripped_fa = os.path.splitext(os.path.basename(gtf_fa))[0] + "-noversions.fa"
    stripped_fa = os.path.join(pizzlydir, stripped_fa)
    gtf_fa = fasta.strip_transcript_versions(gtf_fa, stripped_fa)
    fraglength = get_fragment_length(data)
    cachefile = os.path.join(pizzlydir, "pizzly.cache")
    fusions = kallisto.get_kallisto_fusions(data)
    pizzlypath = config_utils.get_program("pizzly", dd.get_config(data))
    outdir = pizzly(pizzlypath, gtf, gtf_fa, fraglength, cachefile, pizzlydir,
                    fusions, samplename, data)
    return outdir

def pizzly(pizzly_path, gtf, gtf_fa, fraglength, cachefile, pizzlydir, fusions,
           samplename, data):
    outdir = os.path.join(pizzlydir, samplename)
    out_stem = os.path.join(outdir, samplename)
    pizzly_gtf = make_pizzly_gtf(gtf, os.path.join(pizzlydir, "pizzly.gtf"), data)
    sentinel = os.path.join(out_stem, "-flat-filtered.tsv")
    pizzlycalls = out_stem + ".json"
    if not file_exists(pizzlycalls):
        with file_transaction(data, outdir) as tx_out_dir:
            safe_makedir(tx_out_dir)
            tx_out_stem = os.path.join(tx_out_dir, samplename)
            with file_transaction(cachefile) as tx_cache_file:
                cmd = ("{pizzly_path} -k 31 --gtf {pizzly_gtf} --cache {tx_cache_file} "
                    "--align-score 2 --insert-size {fraglength} --fasta {gtf_fa} "
                    "--output {tx_out_stem} {fusions}")
                message = ("Running pizzly on %s." % fusions)
                do.run(cmd.format(**locals()), message)
    flatfile = out_stem + "-flat.tsv"
    filteredfile = out_stem + "-flat-filtered.tsv"
    flatten_pizzly(pizzlycalls, flatfile, data)
    filter_pizzly(flatfile, filteredfile, data)
    return outdir

def make_pizzly_gtf(gtf_file, out_file, data):
    """
    pizzly needs the GTF to be in gene -> transcript -> exon order for each
    gene. it also wants the gene biotype set as the source
    """
    if file_exists(out_file):
        return out_file
    db = gtf.get_gtf_db(gtf_file)
    with file_transaction(data, out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            for gene in db.features_of_type("gene"):
                children = [x for x in db.children(id=gene)]
                for child in children:
                    if child.attributes.get("gene_biotype", None):
                        gene_biotype = child.attributes.get("gene_biotype")
                gene.attributes['gene_biotype'] = gene_biotype
                gene.source = gene_biotype[0]
                print(gene, file=out_handle)
                for child in children:
                    child.source = gene_biotype[0]
                    # gffread produces a version-less FASTA file
                    child.attributes.pop("transcript_version", None)
                    print(child, file=out_handle)
    return out_file

def flatten_pizzly(in_file, out_file, data):
    pizzlyflatten = config_utils.get_program("pizzly_flatten_json.py", data)
    if file_exists(out_file):
        return out_file
    cmd = "{pizzlyflatten} {in_file} > {tx_out_file}"
    message = "Flattening {in_file} to {out_file}."
    with file_transaction(data, out_file) as tx_out_file:
        do.run(cmd.format(**locals()), message.format(**locals()))
    return out_file

def filter_pizzly(in_file, out_file, data):
    df = pd.read_csv(in_file, header=0, sep="\t")
    df = df.query('paircount > 1 and splitcount > 1')
    if file_exists(out_file):
        return out_file
    with file_transaction(out_file) as tx_out_file:
        df.to_csv(tx_out_file, sep="\t", index=False)
    return out_file
