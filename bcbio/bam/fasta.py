from Bio import SeqIO

def sequence_length(fasta):
    """
    return a dict of the lengths of sequences in a fasta file
    """
    file_handle = open(fasta)
    in_handle = SeqIO.parse(file_handle, "fasta")
    records = {record.id: len(record) for record in in_handle}
    file_handle.close()
    return records
