# Detecting gene fusions with bulk RNA-seq data

- `fusion_caller` - Specify a standalone fusion caller for fusion mode. Supports oncofuse for STAR/tophat runs, pizzly and ericscript for all runs. If a standalone caller is specified (i.e. pizzly or ericscript ), fusion detection will not be performed with aligner. oncofuse only supports human genome builds GRCh37 and hg19. ericscript supports human genome builds GRCh37, hg19 and hg38 after installing the associated fusion databases (Customizing data installation).

- `known_fusions` A TAB-delimited file of the format gene1<tab>gene2, where gene1 and gene2 are identifiers of genes specified under gene_name in the attributes part of the GTF file.

