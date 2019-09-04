"""
Wrapper for Salmon:
https://github.com/COMBINE-lab/salmon
http://biorxiv.org/content/early/2015/06/27/021592
"""

import os
import numpy as np

from bcbio.rnaseq import sailfish
import bcbio.pipeline.datadict as dd
from bcbio.utils import (file_exists, safe_makedir, is_gzipped)
import bcbio.utils as utils
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.pipeline import config_utils, fastq
from bcbio import bam
from bcbio.log import logger

def run_salmon_bam(data):
    samplename = dd.get_sample_name(data)
    work_dir = dd.get_work_dir(data)
    salmon_dir = os.path.join(work_dir, "salmon", samplename)
    gtf_file = dd.get_gtf_file(data)
    bam_file = dd.get_transcriptome_bam(data)
    out_file = salmon_quant_bam(bam_file, salmon_dir, gtf_file, data)
    data = dd.set_salmon(data, out_file)
    data = dd.set_salmon_dir(data, salmon_dir)
    data = dd.set_salmon_fraglen_file(data, _get_fraglen_file(salmon_dir))
    data = dd.update_summary_qc(data, "salmon", base=dd.get_salmon_fraglen_file(data))
    return [[data]]

def run_salmon_reads(data):
    data = utils.to_single_data(data)
    files = dd.get_input_sequence_files(data)
    if bam.is_bam(files[0]):
        files = fastq.convert_bam_to_fastq(files[0], data["dirs"]["work"],
                                           data, data["dirs"], data["config"])
    samplename = dd.get_sample_name(data)
    work_dir = dd.get_work_dir(data)
    salmon_dir = os.path.join(work_dir, "salmon", samplename)
    gtf_file = dd.get_gtf_file(data)
    if len(files) == 2:
        fq1, fq2 = files
    else:
        fq1, fq2 = files[0], None
    index = salmon_index(gtf_file, data, os.path.dirname(salmon_dir))
    out_file = salmon_quant_reads(fq1, fq2, salmon_dir, gtf_file, data, index)
    data = dd.set_salmon(data, out_file)
    data = dd.set_salmon_dir(data, salmon_dir)
    data = dd.set_salmon_fraglen_file(data, _get_fraglen_file(salmon_dir))
    data = dd.update_summary_qc(data, "salmon", base=dd.get_salmon_fraglen_file(data))
    return [[data]]

def run_salmon_decoy(data):
    data = utils.to_single_data(data)
    files = dd.get_input_sequence_files(data)
    if bam.is_bam(files[0]):
        files = fastq.convert_bam_to_fastq(files[0], data["dirs"]["work"],
                                           data, data["dirs"], data["config"])
    samplename = dd.get_sample_name(data)
    work_dir = dd.get_work_dir(data)
    salmon_dir = os.path.join(work_dir, "salmon", samplename)
    gtf_file = dd.get_gtf_file(data)
    if len(files) == 2:
        fq1, fq2 = files
    else:
        fq1, fq2 = files[0], None
    index = salmon_decoy_index(gtf_file, data, os.path.dirname(salmon_dir))
    out_file = salmon_quant_reads(fq1, fq2, salmon_dir, gtf_file, data, index)
    data = dd.set_salmon(data, out_file)
    data = dd.set_salmon_dir(data, salmon_dir)
    data = dd.set_salmon_fraglen_file(data, _get_fraglen_file(salmon_dir))
    data = dd.update_summary_qc(data, "salmon", base=dd.get_salmon_fraglen_file(data))
    return [[data]]

def salmon_quant_reads(fq1, fq2, salmon_dir, gtf_file, data, index):
    samplename = dd.get_sample_name(data)
    quant_dir = os.path.join(salmon_dir, "quant")
    safe_makedir(salmon_dir)
    out_file = os.path.join(quant_dir, "quant.sf")
    if file_exists(out_file):
        return out_file
    num_cores = dd.get_num_cores(data)
    strandedness = dd.get_strandedness(data).lower()
    salmon = config_utils.get_program("salmon", dd.get_config(data))
    libtype = sailfish._libtype_string(fq1, fq2, strandedness)
    num_cores = dd.get_num_cores(data)
    resources = config_utils.get_resources("salmon", dd.get_config(data))
    params = ""
    if resources.get("options") is not None:
        params = " ".join([str(x) for x in resources.get("options", [])])
    cmd = ("{salmon} quant {libtype} -i {index} -p {num_cores} "
           " --validateMappings "
           " --seqBias "
           "-o {tx_out_dir} {params} ")
    fq1_cmd = "<(cat {fq1})" if not is_gzipped(fq1) else "<(gzip -cd {fq1})"
    fq1_cmd = fq1_cmd.format(fq1=fq1)
    if not fq2:
        cmd += " -r {fq1_cmd} "
    else:
        cmd += " --gcBias "
        fq2_cmd = "<(cat {fq2})" if not is_gzipped(fq2) else "<(gzip -cd {fq2})"
        fq2_cmd = fq2_cmd.format(fq2=fq2)
        cmd += " -1 {fq1_cmd} -2 {fq2_cmd} "
    # skip --useVBOpt for now, it can cause segfaults
    cmd += "--numBootstraps 30 "
    with file_transaction(data, quant_dir) as tx_out_dir:
        message = ("Quantifying transcripts in %s and %s with Salmon."
                   %(fq1, fq2))
        do.run(cmd.format(**locals()), message, None)
        sailfish.sleuthify_sailfish(tx_out_dir)
    return out_file

def salmon_quant_bam(bam_file, salmon_dir, gtf_file, data):
    samplename = dd.get_sample_name(data)
    quant_dir = os.path.join(salmon_dir, "quant")
    safe_makedir(salmon_dir)
    out_file = os.path.join(quant_dir, "quant.sf")
    if file_exists(out_file):
        return out_file
    if dd.get_transcriptome_fasta(data):
        gtf_fa = dd.get_transcriptome_fasta(data)
    else:
        gtf_fa = sailfish.create_combined_fasta(data)
    num_cores = dd.get_num_cores(data)
    strandedness = dd.get_strandedness(data).lower()
    salmon = config_utils.get_program("salmon", dd.get_config(data))
    libtype = _libtype_string(bam_file, strandedness)
    num_cores = dd.get_num_cores(data)
    cmd = ("{salmon} quant {libtype} -p {num_cores} -t {gtf_fa} "
           "-o {tx_out_dir} -a {bam_file} ")
    cmd += "--numBootstraps 30 "
    with file_transaction(data, quant_dir) as tx_out_dir:
        message = "Quantifying transcripts in %s with Salmon." % bam_file
        do.run(cmd.format(**locals()), message, None)
    return out_file

def _libtype_string(bam_file, strandedness):
    libtype = "-l I" if bam.is_paired(bam_file) else "-l "
    strand = sailfish._sailfish_strand_string(strandedness)
    return libtype + strand

def run_salmon_index(*samples):
    for data in dd.sample_data_iterator(samples):
        work_dir = dd.get_work_dir(data)
        salmon_dir = os.path.join(work_dir, "salmon")
        gtf_file = dd.get_gtf_file(data)
        salmon_index(gtf_file, data, salmon_dir)
    return samples

def salmon_index(gtf_file, data, out_dir):
    out_dir = os.path.join(out_dir, "index", sailfish.get_build_string(data))
    if dd.get_disambiguate(data):
        out_dir = "-".join([out_dir] + dd.get_disambiguate(data))
    salmon = config_utils.get_program("salmon", dd.get_config(data))
    num_cores = dd.get_num_cores(data)
    if dd.get_transcriptome_fasta(data):
        gtf_fa = dd.get_transcriptome_fasta(data)
    else:
        gtf_fa = sailfish.create_combined_fasta(data)
    assert file_exists(gtf_fa), "%s was not found, exiting." % gtf_fa
    tmpdir = dd.get_tmp_dir(data)
    out_file = os.path.join(out_dir, "versionInfo.json")
    if file_exists(out_file):
        logger.info("Transcriptome index for %s detected, skipping building." % gtf_fa)
        return out_dir
    files = dd.get_input_sequence_files(data)
    kmersize = sailfish.pick_kmersize(files[0])
    with file_transaction(data, out_dir) as tx_out_dir:
        cmd = "{salmon} index -k {kmersize} -p {num_cores} -i {tx_out_dir} -t {gtf_fa}"
        message = "Creating Salmon index for {gtf_fa} with {kmersize} bp kmers."
        do.run(cmd.format(**locals()), message.format(**locals()), None)
    return out_dir

def salmon_decoy_index(gtf_file, data, out_dir):
    input_dir = os.path.join(dd.get_work_dir(data), "inputs", "transcriptome")
    decoy_transcriptome = os.path.join(input_dir, sailfish.get_build_string(data) + "-decoy.fa")
    out_dir = os.path.join(out_dir, "index", sailfish.get_build_string(data))
    if dd.get_disambiguate(data):
        out_dir = "-".join([out_dir] + dd.get_disambiguate(data))
    salmon = config_utils.get_program("salmon", dd.get_config(data))
    num_cores = dd.get_num_cores(data)
    if dd.get_transcriptome_fasta(data):
        gtf_fa = dd.get_transcriptome_fasta(data)
    else:
        gtf_fa = sailfish.create_combined_fasta(data)
    assert file_exists(gtf_fa), "%s was not found, exiting." % gtf_fa
    decoy_sequence_file = get_decoy_sequence_file(data)
    decoy_name_file = get_decoy_name_file(data)
    gtf_fa = create_decoy_transcriptome(gtf_fa, get_decoy_sequence_file(data), decoy_transcriptome)
    out_file = os.path.join(out_dir, "versionInfo.json")
    if file_exists(out_file):
        logger.info("Transcriptome index for %s detected, skipping building." % gtf_fa)
        return out_dir
    files = dd.get_input_sequence_files(data)
    kmersize = sailfish.pick_kmersize(files[0])
    with file_transaction(data, out_dir) as tx_out_dir:
        cmd = ("{salmon} index -k {kmersize} -p {num_cores} -i {tx_out_dir} -t {gtf_fa} "
                "--decoys {decoy_name_file} ")
        message = "Creating decoy-aware Salmon index for {gtf_fa} with {kmersize} bp kmers."
        do.run(cmd.format(**locals()), message.format(**locals()), None)
    return out_dir

def create_decoy_transcriptome(transcriptome, decoys, out_file):
    if file_exists(out_file):
        return out_file
    with file_transaction(out_file) as tx_out_file:
        cmd = f"cat {transcriptome} {decoys} > {tx_out_file}"
        message = f"Making decoy transcriptome with {transcriptome} and {decoys} to {out_file}."
        do.run(cmd, message)
    return out_file

def _get_fraglen_file(salmondir):
    flenfile = os.path.join(salmondir, "quant", "libParams", "flenDist.txt")
    if not file_exists(flenfile):
        return None
    else:
        return flenfile

def parse_fragment_length_file(filename):
    with open(filename) as in_handle:
        flens = [float(x) for x in next(in_handle).split("\t")]
    return flens

def estimate_fragment_size(data):
    filename = dd.get_salmon_fraglen_file(data)
    if not file_exists(filename):
        return None
    flen = parse_fragment_length_file(filename)
    return float(np.sum(np.multiply(flen, range(len(flen)))))

def get_decoy_sequence_file(data):
    """
    decoys are located in the rnaseq/salmon-decoys directory
    decoys.fa: decoy sequences to append to transcriptome
    """
    rnaseq_dir = os.path.dirname(dd.get_gtf_file(data))
    decoy_dir = os.path.join(rnaseq_dir, "salmon-decoys")
    decoys = os.path.join(decoy_dir, "decoys.fa")
    if file_exists(decoys):
        return decoys
    else:
        return None

def get_decoy_name_file(data):
    """
    decoys are located in the rnaseq/salmon-decoys directory
    decoys.txt: names of decoy sequences
    """
    rnaseq_dir = os.path.dirname(dd.get_gtf_file(data))
    decoy_dir = os.path.join(rnaseq_dir, "salmon-decoys")
    decoys = os.path.join(decoy_dir, "decoys.txt")
    if file_exists(decoys):
        return decoys
    else:
        return None
