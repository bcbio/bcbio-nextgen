import tempfile
import os
import bcbio.pipeline.datadict as dd
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.utils import file_exists, get_in, safe_makedir, is_gzipped
from bcbio.pipeline import config_utils

def run_sailfish(sample):
    names = sample["rgnames"]
    if len(sample["files"]) == 2:
        fq1, fq2 = sample["files"]
    else:
        fq1, fq2 = sample["files"][0], None
    align_dir = os.path.join(sample["dirs"]["work"], "sailfish", names["sample"])
    safe_makedir(align_dir)
    gtf_file = sample["genome_resources"]["rnaseq"].get("transcripts")
    assert file_exists(gtf_file), "%s was not found, exiting." % gtf_file
    fasta_file = get_in(sample, ("reference", "fasta", "base"))
    assert file_exists(fasta_file), "%s was not found, exiting." % fasta_file
    stranded = get_in(sample["config"], ("algorithm", "strandedness"),
                      "unstranded").lower()
    out_dir = sailfish(fq1, fq2, align_dir, gtf_file, fasta_file, stranded, sample)
    sample["sailfish_dir"] = out_dir
    return [[sample]]


def sailfish(fq1, fq2, align_dir, gtf_file, ref_file, strandedness, data):
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
    return align_dir

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
