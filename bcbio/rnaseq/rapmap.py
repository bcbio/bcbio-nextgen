"""
Wrapper for RapMap:
https://github.com/COMBINE-lab/RapMap
http://biorxiv.org/content/early/2015/10/22/029652
"""
import os

from bcbio.rnaseq import sailfish
import bcbio.pipeline.datadict as dd
from bcbio.utils import (file_exists, safe_makedir, is_gzipped)
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.pipeline import config_utils
from bcbio.ngsalign import postalign
from bcbio import bam
import bcbio.bam.fasta as fasta

def run_rapmap_index(*samples):
    for data in dd.sample_data_iterator(samples):
        work_dir = dd.get_work_dir(data)
        rapmap_dir = os.path.join(work_dir, "rapmap")
        gtf_file = dd.get_gtf_file(data)
        fasta_file = dd.get_ref_file(data)
        assert file_exists(fasta_file), "%s was not found, exiting." % fasta_file
        rapmap_index(gtf_file, fasta_file, "quasi", data, rapmap_dir)
    return samples

def run_rapmap_align(data):
    samplename = dd.get_sample_name(data)
    files = dd.get_input_sequence_files(data)
    work_dir = dd.get_work_dir(data)
    if len(files) == 2:
        fq1, fq2 = files
    else:
        fq1, fq2 = files[0], None
    rapmap_dir = os.path.join(work_dir, "rapmap", samplename)
    gtf_file = dd.get_gtf_file(data)
    fasta_file = dd.get_ref_file(data)
    out_file = rapmap_align(fq1, fq2, rapmap_dir, gtf_file, fasta_file,
                            "quasi", data)
    data = dd.set_transcriptome_bam(data, out_file)
    return [[data]]

def rapmap_index(gtf_file, ref_file, algorithm, data, out_dir):
    out_dir = os.path.join(out_dir, "index")
    valid_indexes = ["pseudoindex", "quasiindex"]
    index_type = algorithm + "index"
    assert index_type in valid_indexes, \
        "RapMap only supports %s indices." % valid_indexes
    out_dir = os.path.join(out_dir, index_type, sailfish.get_build_string(data))
    if dd.get_disambiguate(data):
        out_dir = "-".join([out_dir] + dd.get_disambguate(data))
    rapmap = config_utils.get_program("rapmap", dd.get_config(data))
    # use user supplied transcriptome FASTA file if it exists
    if dd.get_transcriptome_fasta(data):
        gtf_fa = dd.get_transcriptome_fasta(data)
    else:
        gtf_fa = sailfish.create_combined_fasta(data)
    gtf_fa = fasta.strip_transcript_versions(gtf_fa, os.path.join(out_dir, os.path.basename(gtf_fa)))
    tmpdir = dd.get_tmp_dir(data)
    if file_exists(os.path.join(out_dir, "sa.bin")):
        return out_dir
    files = dd.get_input_sequence_files(data)
    kmersize = sailfish.pick_kmersize(files[0])
    message = "Creating rapmap {index_type} for {gtf_fa} with {kmersize} bp kmers."
    with file_transaction(out_dir) as tx_out_dir:
        cmd = "{rapmap} {index_type} -k {kmersize} -i {tx_out_dir} -t {gtf_fa}"
        do.run(cmd.format(**locals()), message.format(**locals()), None)
    return out_dir

def rapmap_align(fq1, fq2, rapmap_dir, gtf_file, ref_file, algorithm, data):
    valid_algorithms = ["pseudo", "quasi"]
    assert algorithm in valid_algorithms, \
        "RapMap algorithm needs to be one of %s." % valid_algorithms
    safe_makedir(rapmap_dir)
    base_dir = os.path.dirname(rapmap_dir)
    samplename = dd.get_sample_name(data)
    out_file = os.path.join(rapmap_dir, samplename + ".bam")
    if file_exists(out_file):
        bam.index(out_file, dd.get_config(data))
        return out_file
    rapmap_index_loc = rapmap_index(gtf_file, ref_file, algorithm, data,
                                    base_dir)
    num_cores = dd.get_num_cores(data)
    algorithm_subcommand = algorithm + "map"
    rapmap = config_utils.get_program("rapmap", dd.get_config(data))
    cmd = "{rapmap} {algorithm_subcommand} -t {num_cores} -i {rapmap_index_loc} "
    fq1_cmd = "{fq1} " if not is_gzipped(fq1) else "<(gzip -cd {fq1}) "
    fq1_cmd = fq1_cmd.format(fq1=fq1)
    if not fq2:
        cmd += "-r {fq1_cmd} "
    else:
        fq2_cmd = "{fq2} " if not is_gzipped(fq2) else "<(gzip -cd {fq2}) "
        fq2_cmd = fq2_cmd.format(fq2=fq2)
        cmd += "-1 {fq2_cmd} -2 {fq2_cmd} "
    with file_transaction(out_file) as tx_out_file:
        cmd += "| " + postalign.sam_to_sortbam_cl(data, tx_out_file)
        run_message = ("%smapping %s and %s to %s with Rapmap. "
                       % (algorithm, fq1, fq2, rapmap_index_loc))
        do.run(cmd.format(**locals()), run_message, None)
    bam.index(out_file, dd.get_config(data))
    return out_file
