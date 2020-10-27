# fast RNA-seq

This mode of `bcbio-nextgen` quantitates transcript expression using [Salmon](https://salmon.readthedocs.io/en/latest/) and does nothing else. It is an order of magnitude faster or more than running the full RNA-seq analysis. The cost of the increased speed is that you will have much less information about your samples at the end of the run, which can make troubleshooting trickier. Invoke with `analysis: fastrna-seq`.


## Parameters
- `transcriptome_fasta` An optional FASTA file of transcriptome sequences
to quantitate rather than using bcbio installed transcriptome sequences.
