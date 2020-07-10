# 3' DGE
3' DGE has some additional complexity when compared to standard bulk RNA-seq experiments. 3' DGE is sequencing just the 3' end of each transcript, so quantification needs to proceed differently. 3' DGE also often incorporates UMIs, which need to be accounted for during the quantification step of RNA-seq. Different
kits have different ways of incorporating the UMI and may or may not include optional well barcodes
for plates and sample barcodes for individual samples. bcbio can be extended to handle arbitrary kits,
and comes with support for several commonly used DGE kits out of the box.

## Description of example dataset
This example takes a small sample of reads generated from the [QIAseq UPX 3' Transcriptome kit](https://www.qiagen.com/us/products/discovery-and-translational-research/next-generation-sequencing/rna-sequencing/three-rnaseq/qiaseq-upx-3-transcriptome-kits/). This particular kit has 96 well and 384 well versions,
which are supported in bcbio via the `umi_type: qiaseq-upx-96` and `umi_type: qiaseq-upx-384` options.

### 1. Download the example data and configuration files
This downloads the input data, creates the project structure and example configuration files.

#### 1.1 Create input directory and download FASTQ files.
```bash
mkdir qiaseq-upx-96-example
cd qiaseq-upx-96-example
mkdir -p fastq
cd fastq
wget --no-check-certificate http://s3.amazonaws.com/bcbio-nextgen/dge_userstory_data/fastq/qiaseq-upx_R1.fastq.gz
wget --no-check-certificate http://s3.amazonaws.com/bcbio-nextgen/dge_userstory_data/fastq/qiaseq-upx_R2.fastq.gz
cd ..
```

#### 1.2 Download template YAML file describing 3' DGE analysis

```bash
wget --no-check-certificate http://s3.amazonaws.com/bcbio-nextgen/dge_userstory_data/qiaseq-upx.yaml
```

qiaseq-upx.yaml: 
```yaml
details:
  - analysis: scrna-seq
    genome_build: hg38
    algorithm:
      umi_type: qiaseq-upx-96
      cellular_barcode_correction: 1
      minimum_barcode_depth: 0
upload:
  dir: ../final
```

#### 1.3 Create a sample sheet

```bash
wget --no-check-certificate http://s3.amazonaws.com/bcbio-nextgen/dge_userstory_data/example_dge.csv
```

example_dge.csv:
```
samplename,description
qiaseq-upx_R1.fastq.gz,testrun
```

### 2. Generate YAML config file for analysis
```
bcbio_nextgen.py -w template qiaseq-upx.yaml example_dge.csv fastq
```

In the result you should see a folder structure:
```
example_dge
|---config
|---final
|---work
```

`example_dge/config/example_dge.yaml` is the main config file to run the bcbio project. You will
see this file has a copy of the parameters in `qiaseq-upx.yaml` for each sample.

### 3. Run the analysis
This will run the analysis on a local machine, using just one core.
```bash
cd example_dge/work
bcbio_nextgen.py ../config/example_dge.yaml -n 1
```

## Parameters

* `umi_type` DGE kit: [harvard-scrb, qiaseq-upx-96, qiaseq-upx-384]
* `minimum_barcode_depth=0` Cellular barcodes with less reads are discarded. This should be 0 if you are using one of these plate-based kits.
* `singlecell_quantifier=rapmap` Quantifier to use for single-cell RNA-sequencing. Supports `rapmap` or `kallisto`. We recommend using `rapmap`.
* *optional* `transcriptome_fasta` alternative transcriptome reference.
* *optional* `transcriptome_gtf` An optional GTF file of the transcriptome to quantitate, rather than the bcbio installed version. 

## Output

Project directory:
```
├── bcbio-nextgen-commands.log -- commands run by bcbio
├── bcbio-nextgen.log -- logging information from bcbio run
├── cb-histogram.txt -- histogram of reads per cellular barcode
├── data_versions.csv -- version information for data used by bcbio
├── metadata.csv -- provided metadata about each sample
├── programs.txt -- program versions of tools run
├── project-summary.yaml -- YAML description of project with derived metadata
├── tagcounts-dupes.mtx -- Matrix Market of gene counts without UMI duplicate removed 
├── tagcounts-dupes.mtx.colnames -- column names to go with tagcounts-dupes.mtx
├── tagcounts-dupes.mtx.rownames -- row names to go with tagcounts-dupes.mtx
├── tagcounts.mtx -- Matrix Market of gene counts, use these for downstream analyses
├── tagcounts.mtx.colnames -- column names to go with tagcounts.mtx
├── tagcounts.mtx.metadata -- optional sample-level metadata for samples in tagcounts.mtx
├── tagcounts.mtx.rownames -- row names to go with tagcounts.mtx
└── transcriptome
    └── mm10.fa -- transcriptome used for quantification

```

Sample directories:
```
testrun/
├── testrun-barcodes-filtered.tsv -- filtered list of cell/well barcodes
├── testrun-barcodes.tsv -- list of cell/well barcodes with number of reads assigned to each barcode
└── testrun-transcriptome.bam -- transcriptome alignments 
```

## Downstream analysis
The starting point for downstream analyses will be the count table of counts per gene per cell/well in the  `tagcounts.mtx`, `tagcounts.mtx.rownames` and
`tagcounts.mtx.colnames` files. You can load these into R using the [readMM function](https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/externalFormats.html) from the [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html) package.

You can use any standard count-based differential RNA-seq differential expression tool to operate on
these count tables such as [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)/[edgeR](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)/[limma](https://bioconductor.org/packages/release/bioc/html/limma.html) and
the analysis will be similar to a bulk RNA-seq experiment. With a large number of samples you will find
making [UMAP plots](https://cran.r-project.org/web/packages/umap/vignettes/umap.html) a useful way to visualize the relationships between your samples. 
