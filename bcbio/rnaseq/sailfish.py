import os
from collections import namedtuple
import pandas as pd

from bcbio import utils
import bcbio.pipeline.datadict as dd
import bcbio.rnaseq.gtf as gtf
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.utils import (file_exists, safe_makedir, is_gzipped,
                         R_package_path, Rscript_cmd)
from bcbio.pipeline import config_utils, disambiguate
from bcbio.bam import fastq
from bcbio import bam

def run_sailfish(data):
    samplename = dd.get_sample_name(data)
    files = dd.get_input_sequence_files(data)
    work_dir = dd.get_work_dir(data)
    if len(files) == 2:
        fq1, fq2 = files
    else:
        fq1, fq2 = files[0], None
    if not fastq.is_fastq(fq1):
        return [[data]]
    sailfish_dir = os.path.join(work_dir, "sailfish", samplename)
    gtf_file = dd.get_gtf_file(data)
    assert file_exists(gtf_file), "%s was not found, exiting." % gtf_file
    fasta_file = dd.get_ref_file(data)
    assert file_exists(fasta_file), "%s was not found, exiting." % fasta_file
    stranded = dd.get_strandedness(data).lower()
    out_file = sailfish(fq1, fq2, sailfish_dir, gtf_file, fasta_file, stranded, data)
    data = dd.set_sailfish(data, out_file)
    data = dd.set_sailfish_dir(data, sailfish_dir)
    return [[data]]

def sailfish(fq1, fq2, sailfish_dir, gtf_file, ref_file, strandedness, data):
    safe_makedir(sailfish_dir)
    samplename = dd.get_sample_name(data)
    quant_dir = os.path.join(sailfish_dir, "quant")
    out_file = os.path.join(quant_dir, "quant.sf")
    if file_exists(out_file):
        return out_file
    build_string = get_build_string(data)
    sailfish_idx = sailfish_index(ref_file, gtf_file, data, build_string)
    num_cores = dd.get_num_cores(data)
    sailfish = config_utils.get_program("sailfish", data["config"])
    cmd = "{sailfish} quant -i {sailfish_idx} -p {num_cores} "
    cmd += _libtype_string(fq1, fq2, strandedness)
    fq1_cmd = "{fq1}" if not is_gzipped(fq1) else "<(gzip -cd {fq1})"
    fq1_cmd = fq1_cmd.format(fq1=fq1)
    if not fq2:
        cmd += " -r {fq1_cmd} "
    else:
        fq2_cmd = "{fq2}" if not is_gzipped(fq2) else "<(gzip -cd {fq2})"
        fq2_cmd = fq2_cmd.format(fq2=fq2)
        cmd += " -1 {fq1_cmd} -2 {fq2_cmd} "
    cmd += "--useVBOpt --numBootstraps 30 "
    cmd += "-o {tx_out_dir}"
    message = "Quantifying transcripts in {fq1} and {fq2}."
    with file_transaction(data, quant_dir) as tx_out_dir:
        do.run(cmd.format(**locals()), message.format(**locals()), None)
        sleuthify_sailfish(tx_out_dir)
    return out_file

def sleuthify_sailfish(sailfish_dir):
    """
    if installed, use wasabi to create abundance.h5 output for use with
    sleuth
    """
    if not R_package_path("wasabi"):
        return None
    else:
        rscript = Rscript_cmd()
        cmd = """{rscript} --vanilla -e 'library("wasabi"); prepare_fish_for_sleuth(c("{sailfish_dir}"))'"""
        do.run(cmd.format(**locals()), "Converting Sailfish to Sleuth format.")
    return os.path.join(sailfish_dir, "abundance.h5")

def create_combined_fasta(data):
    """
    if there are genomes to be disambiguated, create a FASTA file of
    all of the transcripts for all genomes
    """
    out_dir = os.path.join(dd.get_work_dir(data), "inputs", "transcriptome")
    items = disambiguate.split([data])
    fasta_files = []
    for i in items:
        odata = i[0]
        gtf_file = dd.get_gtf_file(odata)
        ref_file = dd.get_ref_file(odata)
        out_file = os.path.join(out_dir, dd.get_genome_build(odata) + ".fa")
        if file_exists(out_file):
            fasta_files.append(out_file)
        else:
            out_file = gtf.gtf_to_fasta(gtf_file, ref_file, out_file=out_file)
            fasta_files.append(out_file)
    out_stem = os.path.join(out_dir, dd.get_genome_build(data))
    if dd.get_disambiguate(data):
        out_stem = "-".join([out_stem] + (dd.get_disambiguate(data) or []))
    combined_file = out_stem + ".fa"
    if file_exists(combined_file):
        return combined_file

    fasta_file_string = " ".join(fasta_files)
    cmd = "cat {fasta_file_string} > {tx_out_file}"
    with file_transaction(data, combined_file) as tx_out_file:
        do.run(cmd.format(**locals()), "Combining transcriptome FASTA files.")
    return combined_file

def create_combined_tx2gene(data):
    out_dir = os.path.join(dd.get_work_dir(data), "inputs", "transcriptome")
    items = disambiguate.split([data])
    tx2gene_files = []
    for i in items:
        odata = i[0]
        gtf_file = dd.get_transcriptome_gtf(odata)
        if not gtf_file:
            gtf_file = dd.get_gtf_file(odata)
        out_file = os.path.join(out_dir, dd.get_genome_build(odata) + "-tx2gene.csv")
        if file_exists(out_file):
            tx2gene_files.append(out_file)
        else:
            out_file = gtf.tx2genefile(gtf_file, out_file, tsv=False)
            tx2gene_files.append(out_file)
    combined_file = os.path.join(out_dir, "tx2gene.csv")
    if file_exists(combined_file):
        return combined_file

    tx2gene_file_string = " ".join(tx2gene_files)
    cmd = "cat {tx2gene_file_string} > {tx_out_file}"
    with file_transaction(data, combined_file) as tx_out_file:
        do.run(cmd.format(**locals()), "Combining tx2gene CSV files.")
    return combined_file

def get_build_string(data):
    build_string = dd.get_genome_build(data)
    if dd.get_disambiguate(data):
        build_string = "-".join([build_string] + (dd.get_disambiguate(data) or []))
    return build_string

def run_sailfish_index(*samples):
    samples = [utils.to_single_data(x) for x in samples]
    Build = namedtuple('Build', ['build', 'ref', 'gtf'])
    builds = {Build(get_build_string(x), dd.get_ref_file(x), dd.get_gtf_file(x))
              for x in samples}
    data = samples[0]
    indexdirs = {}
    for build in builds:
        indexdirs[build.build] = sailfish_index(build.ref, build.gtf, data,
                                                build.build)
    return [[x] for x in samples]

def sailfish_index(gtf_file, ref_file, data, build):
    work_dir = dd.get_work_dir(data)
    out_dir = os.path.join(work_dir, "sailfish", "index", build)
    sailfish = config_utils.get_program("sailfish", data["config"])
    num_cores = dd.get_num_cores(data)
    gtf_fa = create_combined_fasta(data)
    if file_exists(os.path.join(out_dir, "versionInfo.json")):
        return out_dir
    with file_transaction(data, out_dir) as tx_out_dir:
        fq1, _ = dd.get_input_sequence_files(data)
        kmersize = pick_kmersize(fq1)
        cmd = ("{sailfish} index -p {num_cores} -t {gtf_fa} -o {tx_out_dir} "
               "-k {kmersize}")
        message = "Creating sailfish index for {gtf_fa} with {kmersize} bp kmers."
        do.run(cmd.format(**locals()), message.format(**locals()), None)
    return out_dir

def _libtype_string(fq1, fq2, strandedness):
    """
    supports just the Tophat unstranded/firstrand/secondstrand
    """
    libtype = "-l I" if fq2 else "-l "
    strand = _sailfish_strand_string(strandedness)
    return libtype + strand

def _sailfish_strand_string(strandedness):
    return {'unstranded': "U",
            'firststrand': "SR",
            'secondstrand': "SF"}.get(strandedness, "U")

def _sailfish_expression_parser(sailfish_file, samplename):
    col_names = ["name", "length", "effectiveLength", "tpm", "numreads"]
    df = pd.read_csv(sailfish_file, comment="#", header=None, skiprows=1, index_col=0,
                     names=col_names, sep="\t")
    df["sample"] = samplename
    return df

def pick_kmersize(fq):
    """
    pick an appropriate kmer size based off of https://www.biostars.org/p/201474/
    tl;dr version: pick 31 unless the reads are very small, if not then guess
    that readlength / 2 is about right.
    """
    if bam.is_bam(fq):
        readlength = bam.estimate_read_length(fq)
    else:
        readlength = fastq.estimate_read_length(fq)
    halfread = int(round(readlength / 2))
    if halfread >= 31:
        kmersize = 31
    else:
        kmersize = halfread
    if kmersize % 2 == 0:
        kmersize += 1
    return kmersize
