from Bio import SeqIO

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
