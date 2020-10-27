from Bio import SeqIO
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction

def sequence_length(fasta):
    """
    return a dict of the lengths of sequences in a fasta file
    """
    sequences = SeqIO.parse(fasta, "fasta")
    records = {record.id: len(record) for record in sequences}
    return records

def sequence_names(fasta):
    """
    return a list of the sequence IDs in a FASTA file
    """
    sequences = SeqIO.parse(fasta, "fasta")
    records = [record.id for record in sequences]
    return records

def total_sequence_length(fasta):
    """
    return the total length of all sequences in a FASTA file
    """
    return sum([x for x in sequence_length(fasta).values()])

def strip_transcript_versions(fasta, out_file):
    """
    strip transcript versions from a FASTA file. these appear like this:
    >ENST00000434970.2 cdna chromosome:GRCh38:14:22439007:22439015:1 etc
    """
    if file_exists(out_file):
        return out_file
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            with open(fasta) as in_handle:
                for line in in_handle:
                    if line.startswith(">"):
                        out_handle.write(line.split(" ")[0].split(".")[0] + "\n")
                    else:
                        out_handle.write(line)
    return out_file
