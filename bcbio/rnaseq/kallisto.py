"""
kallisto

https://github.com/pachterlab/kallisto
"""
import os
import pandas as pd

import bcbio.pipeline.datadict as dd
from bcbio.rnaseq import sailfish
from bcbio.utils import (file_exists, safe_makedir, chdir)
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.pipeline import config_utils
from bcbio.rnaseq import umi
from bcbio.bam import fasta

def run_kallisto_singlecell(data):
    samplename = dd.get_sample_name(data)
    work_dir = dd.get_work_dir(data)
    kallisto_dir = os.path.join(work_dir, "kallisto", samplename)
    gtf_file = dd.get_gtf_file(data)
    files = dd.get_input_sequence_files(data)
    if len(files) == 2:
        fq1, fq2 = files
    else:
        fq1, fq2 = files[0], None
    assert file_exists(gtf_file), "%s was not found, exiting." % gtf_file
    fasta_file = dd.get_ref_file(data)
    assert file_exists(fasta_file), "%s was not found, exiting." % fasta_file
    out_file = kallisto_singlecell(fq1, kallisto_dir, gtf_file, fasta_file, data)
    data = dd.set_kallisto_quant(data, out_file)
    return [[data]]

def kallisto_singlecell(fq1, kallisto_dir, gtf_file, fasta_file, data):
    samplename = dd.get_sample_name(data)
    quant_dir = os.path.join(kallisto_dir, "quant")
    safe_makedir(kallisto_dir)
    num_cores = dd.get_num_cores(data)
    strandedness = dd.get_strandedness(data).lower()
    kallisto = config_utils.get_program("kallisto", dd.get_config(data))
    # unsure how to estimate from single end data, so go with a reasonable default
    frag_length = 250
    batch_file = umi.convert_to_kallisto(data)
    index = kallisto_index(gtf_file, fasta_file, data, kallisto_dir)
    cmd = ("{kallisto} pseudo --umi "
           "-t {num_cores} -o {tx_out_dir} -b {batch_file} -i {index}")
    with chdir(os.path.dirname(batch_file)):
        with file_transaction(data, quant_dir) as tx_out_dir:
            message = ("Quantifying transcripts with Kallisto.")
            do.run(cmd.format(**locals()), message, None)
    kallisto_table(kallisto_dir, index)
    return quant_dir

def kallisto_index(gtf_file, ref_file, data, out_dir):
    out_dir = os.path.join(out_dir, "index")
    out_stem = dd.get_genome_build(data)
    if dd.get_disambiguate(data):
        out_fname = "-".join([out_fname] + dd.get_disambguate(data))
    index_dir = os.path.join(out_dir, out_stem)
    kallisto = config_utils.get_program("kallisto", dd.get_config(data))
    if dd.get_transcriptome_fasta(data):
        gtf_fa = dd.get_transcriptome_fasta(data)
    else:
        gtf_fa = sailfish.create_combined_fasta(data, index_dir)
    out_file = os.path.join(index_dir, out_stem + ".idx")
    if file_exists(out_file):
        return out_file
    with file_transaction(out_file) as tx_out_file:
        cmd = "{kallisto} index -k 31 -i {tx_out_file} {gtf_fa}"
        message = "Creating Kallisto index for {gtf_fa}."
        do.run(cmd.format(**locals()), message.format(**locals()), None)
    return out_file

def kallisto_table(kallisto_dir, index):
    """
    convert kallisto output to a count table where the rows are
    equivalence classes and the columns are cells
    """
    quant_dir = os.path.join(kallisto_dir, "quant")
    out_file = os.path.join(quant_dir, "matrix.csv")
    if file_exists(out_file):
        return out_file
    tsvfile = os.path.join(quant_dir, "matrix.tsv")
    ecfile = os.path.join(quant_dir, "matrix.ec")
    cellsfile = os.path.join(quant_dir, "matrix.cells")
    fastafile = os.path.splitext(index)[0] + ".fa"
    fasta_names = fasta.sequence_names(fastafile)
    ec_names = get_ec_names(ecfile, fasta_names)
    df = pd.read_table(tsvfile, header=None, names=["ec", "cell", "count"])
    df["ec"] = [ec_names[x] for x in df["ec"]]
    df = df.pivot(index='ec', columns='cell', values='count')
    cellnames = get_cell_names(cellsfile)
    colnames = [cellnames[x] for x in df.columns]
    df.columns = colnames
    df.to_csv(out_file)
    return out_file

def get_ec_names(ecfile, fasta_names):
    """
    convert equivalence classes to their set of transcripts
    """
    df = pd.read_table(ecfile, header=None, names=["ec", "transcripts"])
    transcript_groups = [x.split(",") for x in df["transcripts"]]
    transcripts = []
    for group in transcript_groups:
        transcripts.append(":".join([fasta_names[int(x)] for x in group]))
    return transcripts

def get_cell_names(cellsfile):
    """
    get barcode identifies of cells
    """
    with open(cellsfile) as in_handle:
        return [x.strip() for x in in_handle]
