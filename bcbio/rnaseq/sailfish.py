import tempfile
import os
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.utils import file_exists, get_in, safe_makedir

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
    cmd = "sailfish.sh quant -i {sailfish_idx} "
    cmd += _libtype_string(fq1, fq2, strandedness)
    if not fq2:
        cmd += " -r {fq1} "
    else:
        cmd += " -1 {fq1} -2 {fq2} "
    cmd += " --polya -o {tx_out_dir}"
    message = "Quantifying transcripts in {fq1} and {fq2}."
    with file_transaction(data, align_dir) as tx_out_dir:
        do.run(cmd.format(**locals()), message.format(**locals()), None)
    return align_dir

def sailfish_index(gtf_file, ref_file, data):
    gtf_fa_dirty = _gtf_to_fasta(gtf_file, ref_file, data)
    gtf_fa = _clean_gtf_fa(gtf_fa_dirty, data)
    out_dir = tempfile.mkdtemp(prefix="sailfish_index")
    cmd = "sailfish.sh index -t {gtf_fa} -o {out_dir} -k 25"
    message = "Creating sailfish index for {gtf_fa}."
    do.run(cmd.format(**locals()), message.format(**locals()), None)
    return out_dir

def _libtype_string(fq1, fq2, strandedness):
    type = "PE" if fq2 else "SE"
    # XXX: both orientation and the strand flag could be incorrect for
    # stranded protocols. The Sailfish manual has mode details.
    # this is my best guess for what might work for everything we see
    strand = _sailfish_strand_string(strandedness)
    lstring = "-l \"TYPE={type}:STRAND={strand}"
    if fq2:
        orientation = "><"
        lstring += ":ORIENTATION={orientation}"
    lstring += "\""
    return lstring.format(**locals())


def _sailfish_strand_string(strandedness):
    return {'unstranded': "U",
            'firststrand': "S",
            'secondstrand': "A"}.get(strandedness, "U")


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
