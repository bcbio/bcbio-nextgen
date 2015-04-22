"""
implements support for StringTie, intended to be a drop in replacement for
Cufflinks
http://ccb.jhu.edu/software/stringtie/
http://www.nature.com/nbt/journal/v33/n3/full/nbt.3122.html
manual: http://ccb.jhu.edu/software/stringtie/#contact
"""

import os
import pandas as pd
from bcbio.provenance import do
from bcbio.utils import file_exists, safe_makedir
from bcbio.distributed.transaction import file_transaction
import bcbio.pipeline.datadict as dd

def _stringtie_expression(bam, gtf_file, threads=1, out_dir="."):
    """
    only estimate expression the Stringtie, do not assemble new transcripts
    """
    error_message = "The %s file for %s is missing. StringTie has an error."
    base_cmd = ("stringtie -e -b {out_dir} -p {threads} -G {gtf_file} {bam}")
    transcript_file = os.path.join(out_dir, "t_data.ctab")
    exon_file = os.path.join(out_dir, "e_data.ctab")
    if file_exists(transcript_file):
        return exon_file, transcript_file
    cmd = base_cmd.format(**locals())
    do.run(cmd, "Running Stringtie on %s." % bam)
    assert file_exists(exon_file), error_message % ("exon", exon_file)
    assert file_exists(transcript_file), error_message % ("transcript", transcript_file)
    return exon_file, transcript_file

def run_stringtie_expression(data):
    """
    estimate expression from Stringtie, using the bcbio datadict
    does not do transcriptome assembly
    """
    bam = dd.get_work_bam(data)
    gtf_file = dd.get_gtf_file(data)
    num_cores = dd.get_num_cores(data)
    sample_name = dd.get_sample_name(data)
    out_dir = os.path.join("stringtie", sample_name)
    isoform_fpkm = os.path.join(out_dir, sample_name + ".isoform.fpkm")
    gene_fpkm = os.path.join(out_dir, sample_name + ".fpkm")
    if file_exists(isoform_fpkm) and file_exists(gene_fpkm):
        data = dd.set_cufflinks_dir(data, out_dir)
        data = dd.set_fpkm(data, gene_fpkm)
        data = dd.set_fpkm_isoform(data, isoform_fpkm)
        return data
    with file_transaction(data, out_dir) as tx_out_dir:
        exon_file, transcript_file = _stringtie_expression(bam, gtf_file,
                                                           num_cores, tx_out_dir)
        df = _parse_ballgown(transcript_file)
        _write_fpkms(df, tx_out_dir, sample_name)
    data = dd.set_cufflinks_dir(data, out_dir)
    data = dd.set_fpkm(data, gene_fpkm)
    data = dd.set_fpkm_isoform(data, isoform_fpkm)
    return data

def _write_fpkms(df, out_dir, sample_name):
    transcript_file = os.path.join(out_dir, sample_name + ".isoform.fpkm")
    transcripts = df[["t_name", "FPKM"]]
    transcripts.to_csv(transcript_file, sep="\t", header=False, index=False)
    # sum of transcript FPKM is the gene FPKM
    gene_file = os.path.join(out_dir, sample_name + ".fpkm")
    genes = df[["gene_id", "FPKM"]].groupby(['gene_id']).sum()
    genes.to_csv(gene_file, sep="\t", header=False, index=True)
    return transcript_file, gene_file

def _parse_ballgown(in_file):
   return(pd.DataFrame.from_csv(in_file, header=0, sep="\t"))
