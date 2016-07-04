"""
kallisto

https://github.com/pachterlab/kallisto
"""
import os

from bcbio.rnaseq import sailfish
import bcbio.pipeline.datadict as dd
from bcbio.utils import (file_exists, safe_makedir, chdir)
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.pipeline import config_utils
from bcbio.rnaseq import umi

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
