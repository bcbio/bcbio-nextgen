import tempfile
import os
import shutil
import bcbio.pipeline.datadict as dd
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.utils import (file_exists, safe_makedir, is_gzipped, rbind, partition)
from bcbio.pipeline import config_utils
import pandas as pd

def run_sailfish(data):
    samplename = dd.get_sample_name(data)
    files = dd.get_input_sequence_files(data)
    work_dir = dd.get_work_dir(data)
    if len(files) == 2:
        fq1, fq2 = files
    else:
        fq1, fq2 = files[0], None
    align_dir = os.path.join(work_dir, "sailfish", samplename)
    gtf_file = dd.get_gtf_file(data)
    assert file_exists(gtf_file), "%s was not found, exiting." % gtf_file
    fasta_file = dd.get_ref_file(data)
    assert file_exists(fasta_file), "%s was not found, exiting." % fasta_file
    stranded = dd.get_strandedness(data).lower()
    out_file = sailfish(fq1, fq2, align_dir, gtf_file, fasta_file, stranded, data)
    data = dd.set_sailfish(data, out_file)
    return [[data]]

def sailfish(fq1, fq2, align_dir, gtf_file, ref_file, strandedness, data):
    safe_makedir(align_dir)
    samplename = dd.get_sample_name(data)
    out_file = os.path.join(align_dir, samplename + ".sf")
    if file_exists(out_file):
        return out_file
    sailfish_idx = sailfish_index(gtf_file, ref_file, data)
    num_cores = dd.get_num_cores(data)
    sailfish = config_utils.get_program("sailfish", data["config"])
    cmd = "{sailfish} quant -i {sailfish_idx} -p {num_cores} "
    cmd += _libtype_string(fq1, fq2, strandedness)
    fq1_cmd = "{fq1}" if not is_gzipped(fq1) else "<(gzip -cd {fq1})"
    fq1_cmd = fq1_cmd.format(fq1=fq1)
    if not fq2:
        cmd = " -r {fq1_cmd} "
    else:
        fq2_cmd = "{fq2}" if not is_gzipped(fq2) else "<(gzip -cd {fq2})"
        fq2_cmd = fq2_cmd.format(fq2=fq2)
        cmd += " -1 {fq1_cmd} -2 {fq2_cmd} "
    cmd += "-o {tx_out_dir}"
    message = "Quantifying transcripts in {fq1} and {fq2}."
    with file_transaction(data, align_dir) as tx_out_dir:
        do.run(cmd.format(**locals()), message.format(**locals()), None)
    shutil.move(os.path.join(align_dir, "quant.sf"), out_file)
    return out_file

def sailfish_index(gtf_file, ref_file, data):
    sailfish = config_utils.get_program("sailfish", data["config"])
    gtf_fa_dirty = _gtf_to_fasta(gtf_file, ref_file, data)
    gtf_fa = _clean_gtf_fa(gtf_fa_dirty, data)
    out_dir = tempfile.mkdtemp(prefix="sailfish_index")
    cmd = "{sailfish} index -t {gtf_fa} -o {out_dir} -k 25"
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


def _gtf_to_fasta(gtf_file, ref_file, data):
    gtf_fa = tempfile.NamedTemporaryFile(delete=False, suffix=".fa").name
    with file_transaction(data, gtf_fa) as tx_gtf_fa:
        cmd = "gtf_to_fasta {gtf_file} {ref_file} {tx_gtf_fa}"
        message = "Extracting genomic sequences of {gtf_file}."
        do.run(cmd.format(**locals()), message.format(**locals()), None)
    return gtf_fa

def _clean_gtf_fa(gtf_fa, data):
    """
    convert the gtf_to_fasta sequence names to just the transcript ID
    >1 ENST00000389680 chrM+ 648-1601 -> >ENST00000389680
    """
    out_file = tempfile.NamedTemporaryFile(delete=False, suffix=".fa").name
    with file_transaction(data, out_file) as tx_out_file:
        with open(gtf_fa) as in_handle, open(tx_out_file, "w") as out_handle:
            for line in in_handle:
                if line.startswith(">"):
                    line = ">" + line.split()[1] + "\n"
                out_handle.write(line)
    return out_file

def combine_sailfish(samples):
    work_dir = dd.get_in_samples(samples, dd.get_work_dir)
    dont_combine, to_combine = partition(dd.get_sailfish,
                                         dd.sample_data_iterator(samples), True)
    if not to_combine:
        return samples

    out_file = os.path.join(work_dir, "sailfish", "combined.sf")
    if not file_exists(out_file):
        df = None
        for data in to_combine:
            sailfish_file = dd.get_sailfish(data)
            samplename = dd.get_sample_name(data)
            new_df = _sailfish_expression_parser(sailfish_file, samplename)
            if not df:
                df = new_df
            else:
                df = rbind([df, new_df])
        with file_transaction(out_file) as tx_out_file:
            df.to_csv(tx_out_file, sep="\t", index_label="name")

    updated_samples = []
    for data in dd.sample_data_iterator(samples):
        data = dd.set_sailfish_combined(data, out_file)
        updated_samples.append([data])
    return updated_samples

def _sailfish_expression_parser(sailfish_file, samplename):
    col_names = ["name", "length", "tpm", "numreads"]
    df = pd.io.parsers.read_table(sailfish_file, skiprows=11, header=None,
                                  index_col=0,
                                  names=col_names)
    df["sample"] = samplename
    return df
