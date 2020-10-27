"""Function to counts on the fly the spike in sequences given as a parameter in the yaml file"""
import os
import pandas as pd

# from bcbio.rnaseq import sailfish
import bcbio.pipeline.datadict as dd
from bcbio.utils import (file_exists, safe_makedir, is_gzipped, partition, rbind)
import bcbio.utils as utils
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.pipeline import config_utils
from bcbio.rnaseq import sailfish
from bcbio.bam import fastq
from bcbio.log import logger
# from bcbio import bam

def run_counts_spikein(data):
    return [[counts_spikein(data)]]

def counts_spikein(data):
    data = utils.to_single_data(data)
    samplename = dd.get_sample_name(data)
    work_dir = dd.get_work_dir(data)
    salmon_dir = os.path.join(work_dir, "spikein", samplename)
    fasta_file = dd.get_spikein_fasta(data)
    if not fasta_file:
        return data
    files = dd.get_input_sequence_files(data)
    if len(files) == 2:
        fq1, fq2 = files
    else:
        fq1, fq2 = files[0], None
    assert file_exists(fasta_file), "%s was not found, exiting." % fasta_file
    readlength = fastq.estimate_read_length(fq1)
    if readlength % 2 == 0:
        readlength -= 1
    kmersize = min(readlength, 31)
    logger.info("kmersize used for salmon index at spikein quant: %s" % kmersize)
    kmersize = kmersize if not dd.get_analysis(data).lower() == "smallrna-seq" else 15
    fasta_index = _index_spikein(fasta_file, salmon_dir, data, kmersize)
    out_file = _salmon_quant_reads(fq1, fq2, salmon_dir, fasta_index, data)
    data = dd.set_spikein_counts(data, out_file)
    return data

def _salmon_quant_reads(fq1, fq2, salmon_dir, index, data):
    samplename = dd.get_sample_name(data)
    quant_dir = os.path.join(salmon_dir, "quant")
    safe_makedir(salmon_dir)
    out_file = os.path.join(quant_dir, "quant.sf")
    if file_exists(out_file):
        return out_file
    num_cores = dd.get_num_cores(data)
    salmon = config_utils.get_program("salmon", dd.get_config(data))
    num_cores = dd.get_num_cores(data)
    cmd = ("{salmon} quant --validateMappings -l A -i {index} -p {num_cores} "
           "-o {tx_out_dir} ")
    fq1_cmd = "<(cat {fq1})" if not is_gzipped(fq1) else "<(gzip -cd {fq1})"
    fq1_cmd = fq1_cmd.format(fq1=fq1)
    if not fq2:
        cmd += " -r {fq1_cmd} "
    else:
        fq2_cmd = "<(cat {fq2})" if not is_gzipped(fq2) else "<(gzip -cd {fq2})"
        fq2_cmd = fq2_cmd.format(fq2=fq2)
        cmd += " -1 {fq1_cmd} -2 {fq2_cmd} "
    with file_transaction(data, quant_dir) as tx_out_dir:
        message = ("Quantifying transcripts in %s and %s with Salmon."
                   %(fq1, fq2))
        do.run(cmd.format(**locals()), message, None)
    return out_file

def _index_spikein(fasta, out_dir, data, kmer=31):
    out_dir = safe_makedir(os.path.join(out_dir, "index"))
    salmon = config_utils.get_program("salmon", dd.get_config(data))
    num_cores = dd.get_num_cores(data)
    out_file = os.path.join(out_dir, "versionInfo.json")
    if file_exists(out_file):
        return out_dir
    with file_transaction(out_dir) as tx_out_dir:
        cmd = "{salmon} index -k {kmer} -p {num_cores} -i {tx_out_dir} -t {fasta}"
        message = "Creating Salmon index for {fasta}."
        do.run(cmd.format(**locals()), message.format(**locals()), None)
    return out_dir

def combine_spikein(samples):
    work_dir = dd.get_in_samples(samples, dd.get_work_dir)
    sailfish_dir = os.path.join(work_dir, "spikein")
    dont_combine, to_combine = partition(dd.get_spikein_counts,
                                         dd.sample_data_iterator(samples), True)
    if not to_combine:
        return samples

    tidy_file = os.path.join(sailfish_dir, "spikein.sf")
    if not file_exists(tidy_file):
        logger.info("Combining count files into %s." % tidy_file)
        df = pd.DataFrame()
        for data in to_combine:
            sailfish_file = dd.get_spikein_counts(data)
            samplename = dd.get_sample_name(data)
            new_df = sailfish._sailfish_expression_parser(sailfish_file, samplename)
            if df.empty:
                df = new_df
            else:
                df = rbind([df, new_df])
        df["id"] = df.index
        # some versions of the transcript annotations can have duplicated entries
        df = df.drop_duplicates(["id", "sample"])
        with file_transaction(tidy_file) as tx_out_file:
            df.to_csv(tx_out_file, sep="\t", index_label="name")
        logger.info("Finished combining count files into %s." % tidy_file)

    updated_samples = []
    for data in dd.sample_data_iterator(samples):
        data = dd.set_spikein_counts(data, tidy_file)
        updated_samples.append([data])
    return updated_samples
