# Methylation

Whole genome bisulfite sequencing is supported using
the [bismark2](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) pipeline.
It can be turned on by setting `analysis` to `wgbs-seq`.
Right now we only support the TruSeq Methyl Capture EPIC kit -- if you want to run a different setup,
please [let us know](https://github.com/bcbio/bcbio-nextgen/issues) and we can add support for it.
Consider this a beta feature.

## parameters
- `aligner` supports `bismark`
