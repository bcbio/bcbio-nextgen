import os
import bcbio.pipeline.datadict as dd
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.utils import (file_exists, safe_makedir, is_gzipped, rbind, partition,
                         R_package_path, Rscript_cmd)
from bcbio.pipeline import config_utils, disambiguate
from bcbio.rnaseq import gtf
from bcbio.bam import fastq
from bcbio.log import logger
import pandas as pd
import numpy as np

def run_sailfish(data):
    samplename = dd.get_sample_name(data)
    files = dd.get_input_sequence_files(data)
    work_dir = dd.get_work_dir(data)
    if len(files) == 2:
        fq1, fq2 = files
    else:
        fq1, fq2 = files[0], None
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
    kmer_size = int(fastq.estimate_read_length(fq1))
    if kmer_size < 30:
        # kmer size must be odd
        kmer_size = kmer_size - 5
        kmer_size = kmer_size if kmer_size % 2 else kmer_size - 1
    else:
        kmer_size = 25
    sailfish_idx = sailfish_index(gtf_file, ref_file, data, sailfish_dir, kmer_size)
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
    _sleuthify_sailfish(quant_dir)
    return out_file

def _sleuthify_sailfish(sailfish_dir):
    """
    if installed, use wasabi to create abundance.h5 output for use with
    sleuth
    """
    if not R_package_path("wasabi"):
        return None
    else:
        rscript = Rscript_cmd()
        cmd = """{rscript} -e 'library("wasabi"); prepare_fish_for_sleuth(c("{sailfish_dir}"))'"""
        do.run(cmd.format(**locals()), "Converting Sailfish to Sleuth format.")
    return os.path.join(sailfish_dir, "abundance.h5")

def create_combined_fasta(data, out_dir):
    """
    if there are genomes to be disambiguated, create a FASTA file of
    all of the transcripts for all genomes
    """
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
            out_file = _gtf_to_fasta(gtf_file, ref_file, out_file)
            out_file = _clean_gtf_fa(out_file, out_file)
            fasta_files.append(out_file)
    out_stem = os.path.join(out_dir, dd.get_genome_build(data))
    if dd.get_disambiguate(data):
        out_stem = "-".join([out_stem] + dd.get_disambiguate(data))
    combined_file = out_stem + ".fa"
    if file_exists(combined_file):
        return combined_file

    fasta_file_string = " ".join(fasta_files)
    cmd = "cat {fasta_file_string} > {tx_out_file}"
    with file_transaction(combined_file) as tx_out_file:
        do.run(cmd.format(**locals()), "Combining transcriptome FASTA files.")
    return combined_file

def sailfish_index(gtf_file, ref_file, data, out_dir, kmer_size):
    out_dir = os.path.join(out_dir, "index", dd.get_genome_build(data))
    if dd.get_disambiguate(data):
        out_dir = "-".join([out_dir] + dd.get_disambiguate(data))
    sailfish = config_utils.get_program("sailfish", data["config"])
    num_cores = dd.get_num_cores(data)
    gtf_fa = create_combined_fasta(data, out_dir)
    if file_exists(out_dir + "versionInfo.json"):
        return out_dir
    with file_transaction(out_dir) as tx_out_dir:
        cmd = ("{sailfish} index -p {num_cores} -t {gtf_fa} -o {tx_out_dir} "
               "-k {kmer_size}")
        message = "Creating sailfish index for {gtf_fa}."
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

def _gtf_to_fasta(gtf_file, ref_file, out_file):
    with file_transaction(out_file) as tx_gtf_fa:
        cmd = "gtf_to_fasta {gtf_file} {ref_file} {tx_gtf_fa}"
        message = "Extracting genomic sequences of {gtf_file}."
        do.run(cmd.format(**locals()), message.format(**locals()), None)
    return out_file

def _clean_gtf_fa(gtf_fa, out_file):
    """
    convert the gtf_to_fasta sequence names to just the transcript ID
    >1 ENST00000389680 chrM+ 648-1601 -> >ENST00000389680
    """
    with file_transaction(out_file) as tx_out_file:
        with open(gtf_fa) as in_handle, open(tx_out_file, "w") as out_handle:
            for line in in_handle:
                if line.startswith(">"):
                    line = ">" + line.split()[1] + "\n"
                out_handle.write(line)
    return out_file

def combine_sailfish(samples):
    work_dir = dd.get_in_samples(samples, dd.get_work_dir)
    sailfish_dir = os.path.join(work_dir, "sailfish")
    gtf_file = dd.get_in_samples(samples, dd.get_gtf_file)
    dont_combine, to_combine = partition(dd.get_sailfish,
                                         dd.sample_data_iterator(samples), True)
    if not to_combine:
        return samples

    tidy_file = os.path.join(sailfish_dir, "combined.sf")
    transcript_tpm_file = os.path.join(sailfish_dir, "combined.isoform.sf.tpm")
    gene_tpm_file = os.path.join(sailfish_dir, "combined.gene.sf.tpm")
    tx2gene = os.path.join(sailfish_dir, "tx2gene.csv")
    if not all([file_exists(x) for x in [gene_tpm_file, tidy_file,
                                         transcript_tpm_file, tx2gene]]):
        logger.info("Combining count files into %s." % tidy_file)
        df = pd.DataFrame()
        for data in to_combine:
            sailfish_file = dd.get_sailfish(data)
            samplename = dd.get_sample_name(data)
            new_df = _sailfish_expression_parser(sailfish_file, samplename)
            if df.empty:
                df = new_df
            else:
                df = rbind([df, new_df])
        df["id"] = df.index
        # some versions of the transcript annotations can have duplicated entries
        df = df.drop_duplicates(["id", "sample"])
        with file_transaction(tidy_file) as tx_out_file:
            df.to_csv(tx_out_file, sep="\t", index_label="name")
        with file_transaction(transcript_tpm_file) as  tx_out_file:
            df.pivot("id", "sample", "tpm").to_csv(tx_out_file, sep="\t")
        with file_transaction(gene_tpm_file) as  tx_out_file:
            pivot = df.pivot("id", "sample", "tpm")
            tdf = pd.DataFrame.from_dict(gtf.transcript_to_gene(gtf_file),
                                         orient="index")
            tdf.columns = ["gene_id"]
            pivot = pivot.join(tdf)
            pivot = pivot.groupby("gene_id").agg(np.sum)
            pivot.to_csv(tx_out_file, sep="\t")
        tx2gene = gtf.tx2genefile(gtf_file, tx2gene)
        logger.info("Finished combining count files into %s." % tidy_file)

    updated_samples = []
    for data in dd.sample_data_iterator(samples):
        data = dd.set_sailfish_tidy(data, tidy_file)
        data = dd.set_sailfish_transcript_tpm(data, transcript_tpm_file)
        data = dd.set_sailfish_gene_tpm(data, gene_tpm_file)
        data = dd.set_tx2gene(data, tx2gene)
        updated_samples.append([data])
    return updated_samples

def _sailfish_expression_parser(sailfish_file, samplename):
    col_names = ["name", "length", "effectiveLength", "tpm", "numreads"]
    df = pd.read_csv(sailfish_file, comment="#", header=None, skiprows=1, index_col=0,
                     names=col_names, sep="\t")
    df["sample"] = samplename
    return df
