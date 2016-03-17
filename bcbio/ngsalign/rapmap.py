import os
import bcbio.pipeline.datadict as dd
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.utils import (file_exists, safe_makedir, is_gzipped)
from bcbio.pipeline import config_utils
from bcbio.rnaseq.sailfish import create_combined_fasta

def run_rapmap_pseudoalign(data):
    samplename = dd.get_sample_name(data)
    files = dd.get_input_sequence_files(data)
    work_dir = dd.get_work_dir(data)
    if len(files) == 2:
        fq1, fq2 = files
    else:
        fq1, fq2 = files[0], None
    rapmap_dir = os.path.join(work_dir, "rapmap", samplename)
    gtf_file = dd.get_gtf_file(data)
    assert file_exists(gtf_file), "%s was not found, exiting." % gtf_file
    fasta_file = dd.get_ref_file(data)
    assert file_exists(fasta_file), "%s was not found, exiting." % fasta_file
    out_file = rapmap_pseudoalign(fq1, fq2, rapmap_dir, gtf_file, fasta_file,
                                  data)
    data = dd.set_work_bam(data, out_file)
    data = dd.set_transcriptome_bam(data, out_file)
    return data

def rapmap_pseudoalign(fq1, fq2, rapmap_dir, gtf_file, ref_file, data):
    safe_makedir(rapmap_dir)
    samplename = dd.get_sample_name(data)
    out_file = os.path.join(rapmap_dir, samplename + ".bam")
    if file_exists(out_file):
        return out_file
    rapmap_idx = rapmap_index(gtf_file, ref_file, data, sailfish_dir)
    num_cores = dd.get_num_cores(data)
    rapmap = config_utils.get_program("rapmap", data["config"])
    cmd = ("{rapmap} pseudomap -i {rapmap_idx} -t {num_cores} ")
    fq1_cmd = "{fq1}" if not is_gzipped(fq1) else "<(gzip -cd {fq1})"
    fq1_cmd = fq1_cmd.format(fq1=fq1)
    if not fq2:
        cmd += " -r {fq1_cmd} "
    else:
        fq2_cmd = "{fq2}" if not is_gzipped(fq2) else "<(gzip -cd {fq2})"
        fq2_cmd = fq2_cmd.format(fq2=fq2)
        cmd += " -1 {fq1_cmd} -2 {fq2_cmd} "
    message = "pseudomapping transcripts in {fq1} and {fq2}."
    with file_transaction(data, out_file) as tx_out_file:
        cmd += " | " + postalign.sam_to_sortbam_cl(data, tx_out_file)
        do.run(cmd.format(**locals()), message.format(**locals()), None)
    return out_file

def rapmap_index(gtf_file, ref_file, data, out_dir):
    out_dir = os.path.join(out_dir, "index", dd.get_genome_build(data))
    if dd.get_disambiguate(data):
        out_dir = "-".join([out_dir] + dd.get_disambiguate(data))
    rapmap = config_utils.get_program("rapmap", data["config"])
    gtf_fa = create_combined_fasta(data, out_dir)
    if file_exists(out_dir + "rapidx.jfhash"):
        return out_dir
    with file_transaction(out_dir) as tx_out_dir:
        cmd = "{rapmap} pseudoindex -i {tx_out_dir} -t {gtf_fa}"
        message = "Creating RapMap pseudoindex for {gtf_fa}."
        do.run(cmd.format(**locals()), message.format(**locals()), None)
    return out_dir
