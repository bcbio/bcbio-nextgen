# Counting cells and transcripts for inDrops3 data

## Workflow

Bcbio installation paths in this workflow correspond to [O2 bcbio installation](https://wiki.rc.hms.harvard.edu/display/O2).
Adjust to bcbio installation you are working with.

### 1. Check reference genome and transcriptome - is it a mouse project?
- mm10 reference genome: /n/shared_db/bcbio/biodata/genomes/Mmusculus/mm10
- transcriptome_fasta: /n/shared_db/bcbio/biodata/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.fa
- transcriptome_gtf: /n/shared_db/bcbio/biodata/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.gtf

### 2. Create bcbio project structure in /scratch
```
mkdir sc_mouse
cd sc_mouse
mkdir config input final work
```

### 3. Prepare fastq input in sc_mouse/input
- some FC come in 1..4 lanes, merge lanes for every read:
```
cat lane1_r1.fq.gz lane2_r1.fq.gz > project_1.fq.gz
cat lane1_r2.fq.gz lane2_r2.fq.gz > project_2.fq.gz
```
- cat'ing gzip files sounds ridiculous, but works for the most part, for purists:
```
zcat KM_lane1_R1.fastq KM_lane2_R1.fastq.gz | gzip > KM_1.fq.gz
```

- some cores send bz2 files not gz
```
bunzip2 *.bz2
cat *R1.fastq | gzip > sample_1.fq.gz
```

- some cores produce R1,R2,R3,R4, others R1,R2,I1,I2, rename them
```
bcbio_R1 = R1 = 86 or 64 bp transcript read
bcbio_R2 = I1 = 8 bp part 1 of cell barcode
bcbio_R3 = I2 = 8 bp sample (library) barcode
bcbio_R4 = R2 = 14 bp = 8 bp part 2 of cell barcode + 6 bp of transcript UMI
```
- files in sc_mouse/input should be (KM here is project name):
```
KM_1.fq.gz
KM_2.fq.gz
KM_3.fq.gz
KM_4.fq.gz
```

### 4. Specify sample barcodes
Sample barcodes should be in *sc_mouse/config/sample_barcodes.csv*.    

Check out if the sample barcodes provided match the actual barcodes in the data.

```shell
gunzip -c FC_X_3.fq.gz | awk '{if(NR%4 == 2) print $0}' | head -n 400000 | sort | uniq -c | sort -k1,1rn | awk '{print $2","$1}' | head

AGGCTTAG,112303
ATTAGACG,95212
TACTCCTT,94906
CGGAGAGA,62461
CGGAGATA,1116
CGGATAGA,944
GGGGGGGG,852
ATTAGACC,848
ATTAGCCG,840
ATTATACG,699
```

Sometimes you need to reverse complement sample barcodes:
```
cat barcodes_original.csv | awk -F ',' '{print $1}' | tr ACGTacgt TGCAtgca | rev
```

sample_barcodes.csv
```
TCTCTCCG,S01
GCGTAAGA,S02
CCTAGAGT,S03
TCGACTAG,S04
TTCTAGAG,S05
```

### 5. Create bcbio yaml config file
*sc_mouse/config/sc-mouse.yaml*:

```
details:
- algorithm:
    cellular_barcode_correction: 1
    minimum_barcode_depth: 1000
    sample_barcodes: /full/path/sc_mouse/config/sample_barcodes.csv
    transcriptome_fasta: /n/shared_db/bcbio/biodata/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.fa
    transcriptome_gtf: /n/shared_db/bcbio/biodata/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.gtf
    umi_type: harvard-indrop-v3
  analysis: scRNA-seq
  description: PI_name
  files:
  - /full/path/sc_mouse/input/KM_1.fq.gz
  - /full/path/sc_mouse/input/KM_2.fq.gz
  - /full/path/sc_mouse/input/KM_3.fq.gz
  - /full/path/sc_mouse/input/KM_4.fq.gz
  genome_build: mm10
  metadata: {}
fc_name: sc-mouse
upload:
  dir: /full/path/sc_mouse/final
```
Use `cd sc_mouse/input; readlink -f *` to grab full path to each file and paste into yaml.

### 6. Create a batch script
*sc_mouse/config/bcbio.sh*:

```shell
#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=priority        # Partition (queue)
#SBATCH --time=10-00:00             # Runtime in D-HH:MM format
#SBATCH --job-name=km            # Job name
#SBATCH -c 20
#SBATCH --mem-per-cpu=5G            # Memory needed per CPU
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

bcbio_nextgen.py ../config/sc-mouse.yaml -n 20
```
- most projects take < 5days, but some large 4 lane could take more, like 7-8

### 7. Run bcbio

```shell
cd sc_mouse_work
sbatch ../config/bcbio.sh
```

### 8.1 Create Seurat object
Reading matrices from bcbio and creating seurat objects might take a lot of RAM.
If you can't fit your dataset on a laptop, create R conda environment on a cluster,
as [described](https://github.com/hbc/knowledgebase/blob/master/scrnaseq/Single-Cell-conda.md)

```
# O2: use interactive job with 20G RAM
# sometimes Rscript is not working but works from R
# conda activate r
# which R
# Rscript 00_create_seurat_object.R
# conda deactivate
# should have 3 files from bcbio
library(R.utils)

file.rename("tagcounts.mtx", "matrix.mtx")
file.rename("tagcounts.mtx.rownames", "features.tsv")
file.rename("tagcounts.mtx.colnames", "barcodes.tsv")

gzip("matrix.mtx")
gzip("features.tsv")
gzip("barcodes.tsv")

library(Seurat)
counts <- Read10X(data.dir = ".", gene.column = 1)
seurat_object <- CreateSeuratObject(counts = counts, min.features = 100)
saveRDS(seurat_object, "seurat.bcbio.RDS")
```

### 8.2. Create SingleCellExperiment object
You may convert Seurat object to SingleCellExperiment, but also you may choose to save duplicated matrices in SCE.

```
# run:
# interactive session with 20G of RAM
# cd to project/final/project
# conda activate r
# Rscript 01.bcbio2sse.R

library(SingleCellExperiment)
library(Matrix)
library(AnnotationHub)
library(tidyverse)

species = "Mus musculus"
counts = readMM(file.path("tagcounts.mtx"))
dupcounts = readMM(file.path("tagcounts-dupes.mtx"))
rownames = read.csv(file.path("tagcounts.mtx.rownames"), header = F)[["V1"]]
rownames = as.character(rownames)
colnames = read.csv(file.path("tagcounts.mtx.colnames"), header = F)[["V1"]]
colnames = make.names(as.character(colnames))
reads = read.csv(file.path("cb-histogram.txt"), header = F, sep="\t", row.names = 1)
rownames(reads) = make.names(rownames(reads))

counts =  as(counts, "dgCMatrix")
rownames(counts) = rownames
colnames(counts) = colnames
metadata = read.csv(file.path("tagcounts.mtx.metadata"))
rownames(metadata) = colnames
metadata[["nUMI"]] = colSums(counts)
metadata[["nGenes"]] = colSums(counts>0)
metadata[["log10GenesPerUMI"]] = log10(metadata$nGene) / log10(metadata$nUMI)
metadata[["nReads"]] = reads[colnames,]
metadata[["saturation_rate"]] = 1-(colSums(counts)/colSums(dupcounts))
metadata[["dupReads"]] = colSums(dupcounts)
metadata[["dupMeanReads"]] = colMeans(dupcounts)

# check if file is empty and skip if the case
# annotation was download from ensembl biomart to match the version GRCh38.92
# AnnotationHub can be used.

## Load the annotation resource.
ah <- AnnotationHub()
ahDb <- query(ah, pattern= c(species, "EnsDb") )
ahEdb <- ahDb[[rev(names(ahDb))[1]]] # last one is chosen
rows = genes(ahEdb) %>%
    as.data.frame() %>%
    janitor::clean_names() %>%
    dplyr::select(gene_id,
                  gene_name,
                  description,
                  biotype = gene_biotype,
                  entrezid,
                  chrom = seqnames) %>%
    group_by(gene_id, biotype, chrom, description) %>%
    summarise(gene_name = paste(unique(gene_name), collapse = ","),
              entrezid = paste(unique(entrezid), collapse = ",")) %>%
    mutate(gene_name=ifelse(gene_name=="", gene_id, gene_name)) %>%
    as.data.frame()

# mit
rrna = rows %>% dplyr::filter(chrom == "MT") %>% .[["gene_id"]] %>% intersect(., rownames)

metadata[["mtUMI"]] = colSums(counts[rrna,], na.rm = T)
metadata[["mtUMI"]][is.na(metadata[["mtUMI"]])] = 0
metadata[["mitoRatio"]] = metadata$mtUMI/metadata$nUMI

se = SingleCellExperiment(assays=list(raw = counts), colData = metadata)
saveRDS(se, "se.RDS")
```

### 1a. (Optional).
If you care, download fresh transcriptome annotation from Gencode (https://www.gencodegenes.org/mouse/)
(it has chrom names with chr matching mm10 assembly).
```
cd sc_mouse/input
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz
gunzip gencode.vM23.annotation.gtf.gz
gffread -g /n/shared_db/bcbio/biodata/genomes/Mmusculus/mm10/seq/mm10.fa gencode.vM23.annotation.gtf -x gencode.vM23.annotation.cds.fa
```
update sc_mouse/config/sc_mouse.yaml:
```
transcriptome_fasta: gencode.vM23.annotation.cds.fa
transcriptome_gtf: gencode.vM23.annotation.gtf
```

## Parameters

* `umi_type` Single cell library type: [harvard-indrop, harvard-indrop-v2, 10x_v2, icell8, surecell].
* `minimum_barcode_depth=10000` Cellular barcodes with less reads are discarded.
* `sample_barcodes` A file with one sample barcode per line.
If the file contains sample name for each barcode, this will be used to create a `tagcounts.mtx.metadata`
that match each cell with the sample name associated with the barcode.
For inDrops3 protocol sample barcodes are in the fastq file for read3. Example of `sample_barcodes` file:
    ```
    AATTCCGG,sample1
    CCTTGGAA,sample2
    ```
* `singlecell_quantifier=rapmap` Quantifier to use for single-cell RNA-sequencing. Supports `rapmap` or `kallisto`.
* *optional* `cellular_barcodes` A file or a list [file1.txt, file2.txt] of valid cellular barcodes.
* *optional* `cellular_barcode_correction=1` Number of errors to correct in identified cellular barcodes. Requires `cellular_barcodes` option set. Set to 0 to turn off error correction.
* *optional* `transcriptome_fasta` alternative transcriptome reference.
* *optional* `transcriptome_gtf` An optional GTF file of the transcriptome to quantitate, rather than the bcbio installed version. This is recommended for single-cell RNA-sequencing experiments.
* `demultiplexed` If set to True, each file will be treated as a cell or well and not a collection of cells. Use this if your data has already been broken up into cells or wells.

## Output

Project directory `bcbio_run/final/project`:
* `tagcounts.mtx` -- count matrix compatible with dgCMatrix type in R.
* `tagcounts-dupes.mtx` -- count matrix compatible with dgCMatrix type in R but with the duplicated reads counted.
* `tagcounts.mtx.colnames` -- cell names that would be the columns for the matrix.
* `tagcounts.mtx.rownames` -- gene names that would be the rows for the matrix.
* `tagcounts.mtx.metadata` -- metadata that match the colnames for the matrix. This is coming from the barcode.csv file and the metadata given in the YAML config file. for the matrix.
* `cb-histogram.txt` -- total number of dedup reads assigned to a cell. Comparing colSums(tagcounts.mtx) to this number can tell you how many reads mapped to genes.

Sample directories `bcbio_run/final/sample`:
* `SAMPLE-transcriptome.bam` -- BAM file aligned to transcriptome.
* `SAMPLE-mtx.*` -- gene counts as explained in the project directory.

## Steps

Minor steps and exact file locations are omitted.

* = repeated for every sample in sample_barcodes.csv

*bcbio.yaml config*:

```yaml
details:
- algorithm:
    cellular_barcode_correction: 1
    minimum_barcode_depth: 1000
    sample_barcodes: /path/project/config/barcodes.csv
    transcriptome_fasta: /path/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.fa
    transcriptome_gtf: /path/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.gtf
    umi_type: harvard-indrop-v3
  analysis: scRNA-seq
  description: project
  files:
  - /path/project_1.fq.gz
  - /path/project_2.fq.gz
  - /path/project_3.fq.gz
  - /path/project_4.fq.gz
  genome_build: mm10
fc_name: sc-mouse
upload:
  dir: /path/project/final
```

1. parse barcode information from reads 2,3,4 to the fastq read name, [CELL]\_[value]:[UMI]\_[value]:[sample]\_[value]. umis supports many protocols, however, the downside is speed - this step can take up to 3-4 days:
    ```shell
    python umis fastqtransform \
    --separate_cb umis/harvard-indrop-v3-transform.json --cores 16  \
    project_1.fq.gz project_2.fq.gz project_3.fq.gz project_4.fq.gz | \
    seqtk seq -L 20 - | gzip > project.umitransformed.fq.gz
    ```
2. create a fastq file for each sample:
    ```shell
    python umis demultiplex_samples \
    --nedit 1 --barcodes sample_barcodes.csv \
    --out_dir demultiplexed project.umitransformed.fq.gz
    ```
3. Cellular barcode filter
    ```shell
    python umis cb_filter \
    --cores 16 --bc1 harvard-indrop-v3-cb1.txt.gz --nedit 1 \
    --bc2 harvard-indrop-v3-cb2.txt.gz demultiplexed/[sample-barcodeAATTTTT].fq  | \
    gzip -c > project-sample_barcode.filtered.fq.gz
    ```
4. \* create cellular barcode histogram, also creates cb-histogram-filtered.txt for cells with nreads > minimum_barcode_depth:
    ```shell
    python umis cb_histogram project-sample_barcode.filtered.fq.gz > cb-histogram.txt
    ```
5. \* create index genome for rapmap:
    ```shell
    rapmap quasiindex -k 31 -i mm10 -t mm10/ref-transcripts.fa
    ```
6. \* align reads with rapmap:
    ```shell
    rapmap quasimap -t 16 -i mm10 \
    -r <(gzip -cd project-sample_barcode.filtered.fq.gz) | \
    samtools sort -@ 16 -m 1G  -T project-sample_barcode-sorttmp \
    -o project-sample_barcode.bam /dev/stdin
    samtools index -@ 16 project-sample_barcode.bam project-sample_barcode.bam.bai
    ```
7. \* count transcripts:
    ```shell
    python umis fasttagcount --cb_cutoff 1000 \
    --genemap ref-transcripts-tx2gene.tsv
    --cb_histogram project-sample_barcode/cb-histogram.txt \
    --umi_matrix project-sample_barcode-dupes.mtx.full \
    project-sample_barcode.bam project-sample_barcode.mtx.full

    python umis sparse project-sample_barcode.mtx.full \
    project-sample_barcode.mtx

    python umis sparse project-sample_barcode-dupes.mtx.full \
    project-sample_barcode-dupes.mtx
    ```
8. Concatenate all cb-histogram-filtered.txt files:
    ```shell
    cat project-[all-barcodes]/cb-histogram-filtered.txt > cb-histogram.txt
    ```
## Description

bcbio-nextgen supports universal molecular identifiers (UMI) based single-cell RNA-seq analyses. If your single-cell prep does not use universal molecular identifiers (UMI), you can most likely just run the standard RNA-seq pipeline and use the results from that. The UMI are used to discard reads which are possibly PCR duplicates and is very helpful for removing some of the PCR duplicate noise that can dominate single-cell experiments.

Unlike the standard RNA-seq pipeline, the single-cell pipeline expects the FASTQ input files to not be separated by cellular barcode, so each file is a mix of cells identified by a cellular barcode (CB), and unique reads from a transcript are identified with a UMI. bcbio-nextgen inspects each read, identifies the cellular barcode and UMI and puts them in the read name. Then the reads are aligned to the transcriptome with [RapMap](https://github.com/COMBINE-lab/RapMap) and the number of reads aligning to each transcript is counted for each cellular barcode. The output is a table of counts with transcripts as the rows and columns as the cellular barcodes for each input FASTQ file.

Optionally the reads can be quantitated with `kallisto` to output transcript compatibility counts rather than counts per gene ([TCC paper](https://doi.org/10.1186/s13059-016-0970-8)).

To extract the UMI and cellular barcodes from the read, bcbio-nextgen needs to know where the UMI and the cellular barcode are expected to be in the read. Currently there is support for two schemes, the inDrop system from the Harvard single-cell core facility and CEL-seq. If bcbio-nextgen does not support your UMI and barcoding scheme, please open up an issue and we will help implement support for it.

Most of the heavy lifting for this part of bcbio-nextgen is implemented in the [umis](https://github.com/vals/umis) repository.

## References
- [Indrops3 library structure](https://singlecellcore.hms.harvard.edu/resources)
- [Even shorter guide](https://github.com/bcbio/bcbio-nextgen/blob/master/config/templates/indrop-singlecell.yaml)
- [Much more comprehensive guide](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/01_bcbio_run.md)
