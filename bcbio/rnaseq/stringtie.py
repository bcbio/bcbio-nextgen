"""
implements support for StringTie, intended to be a drop in replacement for
Cufflinks
http://ccb.jhu.edu/software/stringtie/
http://www.nature.com/nbt/journal/v33/n3/full/nbt.3122.html
manual: http://ccb.jhu.edu/software/stringtie/#contact
"""

import os
import pandas as pd
import subprocess
import contextlib
from distutils.version import LooseVersion

from bcbio.provenance import do
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
import bcbio.pipeline.datadict as dd

def _stringtie_expression(bam, data, out_dir="."):
    """
    only estimate expression the Stringtie, do not assemble new transcripts
    """
    gtf_file = dd.get_gtf_file(data)
    num_cores = dd.get_num_cores(data)
    error_message = "The %s file for %s is missing. StringTie has an error."
    stringtie = config_utils.get_program("stringtie", data, default="stringtie")
    # don't assemble transcripts unless asked
    exp_flag = ("-e" if "stringtie" not in dd.get_transcript_assembler(data)
                else "")
    base_cmd = ("{stringtie} {exp_flag} -b {out_dir} -p {num_cores} -G {gtf_file} "
                "-o {out_gtf} {bam}")
    transcript_file = os.path.join(out_dir, "t_data.ctab")
    exon_file = os.path.join(out_dir, "e_data.ctab")
    out_gtf = os.path.join(out_dir, "stringtie-assembly.gtf")
    if file_exists(transcript_file):
        return exon_file, transcript_file, out_gtf
    cmd = base_cmd.format(**locals())
    do.run(cmd, "Running Stringtie on %s." % bam)
    assert file_exists(exon_file), error_message % ("exon", exon_file)
    assert file_exists(transcript_file), error_message % ("transcript", transcript_file)
    return transcript_file

def run_stringtie_expression(data):
    """
    estimate expression from Stringtie, using the bcbio datadict
    does not do transcriptome assembly
    """
    bam = dd.get_work_bam(data)
    sample_name = dd.get_sample_name(data)
    out_dir = os.path.join("stringtie", sample_name)
    isoform_fpkm = os.path.join(out_dir, sample_name + ".isoform.fpkm")
    gene_fpkm = os.path.join(out_dir, sample_name + ".fpkm")
    assembly = os.path.abspath(os.path.join(out_dir, "stringtie-assembly.gtf"))
    if file_exists(isoform_fpkm) and file_exists(gene_fpkm):
        data = dd.set_stringtie_dir(data, out_dir)
        data = dd.set_fpkm(data, gene_fpkm)
        data = dd.set_fpkm_isoform(data, isoform_fpkm)
        if "stringtie" in dd.get_transcript_assembler(data):
            assembled_gtfs = dd.get_assembled_gtf(data)
            assembled_gtfs.append(assembly)
            data = dd.set_assembled_gtf(data, assembled_gtfs)
        return data
    with file_transaction(data, out_dir) as tx_out_dir:
        transcript_file = _stringtie_expression(bam, data, tx_out_dir)
        df = _parse_ballgown(transcript_file)
        _write_fpkms(df, tx_out_dir, sample_name)
    data = dd.set_stringtie_dir(data, out_dir)
    data = dd.set_fpkm(data, gene_fpkm)
    data = dd.set_fpkm_isoform(data, isoform_fpkm)
    if "stringtie" in dd.get_transcript_assembler(data):
        assembled_gtfs = dd.get_assembled_gtf(data)
        assembled_gtfs.append(assembly)
        data = dd.set_assembled_gtf(data, assembled_gtfs)
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
    return(pd.read_csv(in_file, header=0, sep="\t"))

def merge(to_merge, ref_file, gtf_file, num_cores, data):
    stringtie = config_utils.get_program("stringtie", data, default="stringtie")
    gtf_list = " ".join(to_merge)
    out_dir = os.path.join("assembly", "stringtie-merge")
    merged_file = os.path.join(out_dir, "merged.gtf")
    cmd = "{stringtie} --merge -o {tx_merged_file} -G {gtf_file} {gtf_list}"
    if not file_exists(merged_file):
        with file_transaction(merged_file) as tx_merged_file:
            message = "Merging transcriptome assemblies with Stringtie."
            do.run(cmd.format(**locals()), message)
    return merged_file

def version(data):
    stringtie = config_utils.get_program("stringtie", data, default="stringtie")
    cmd = "%s --version" % stringtie
    subp = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            shell=True)
    with contextlib.closing(subp.stdout) as stdout:
        for line in stdout:
            version = line.decode().strip()
    return LooseVersion(version)

def supports_merge(data):
    """
    1.2.0 and up supports the --merge option and a planned utility gffcompare
    (https://github.com/gpertea/stringtie/issues/29) will add the class codes
    which will remove the need for cufflinks merge
    """
    gffcompare_installed = config_utils.program_installed("gffcompare", data)
    return version(data) >= LooseVersion("1.2.0") and gffcompare_installed
