## Outputs

bcbio-nextgen runs in a temporary work directory which contains a number of processing intermediates. Pipeline completion extracts the final useful output files into a separate directory, specified by the [Upload](contents/configuration:upload). This configuration allows upload to local directories, Galaxy, or Amazon S3. Once extracting and confirming the output files, you can delete the temporary directory to save space.

### Common files

The output directory contains sample specific output files labeled by sample name and a more general project directory. The sample directories contain all of the sample specific output files, while the project directory contains global files like project summaries or batched population level variant calls. See the [Teaching](teaching) documentation for a full variant calling example with additional details about configuration setting and resulting output files.

Project directory:
* `project-summary.yaml` -- Top level YAML format summary file with statistics on read alignments and duplications as well as analysis specific metrics.
* `programs.txt` -- Program versions for bcbio-nextgen and software run in the pipeline. This enables reproduction of analyses.
* `multiqc` run [MultiQC](https://multiqc.info/) to gather all QC metrics from different tools, such as, cutadapt, featureCounts, samtools, STAR .. into an unique HTML report.
* `metadata.csv` -- CSV with the metadata in the YAML file.
* `data_versions.csv` -- Data versions for bcbio-nextgen and software

Sample directories:
* `SAMPLE/qc` -- Directory of quality control runs for the sample. These include charts and metrics for assessing quality of sequencing and analysis.
* `SAMPLE-ready.bam` -- A prepared BAM file of the aligned reads. Depending on the analysis used, this may include trimmed, recalibrated and realigned reads following alignment.

### Variant calling

Project directory:
* `grading-summary.csv` -- Grading details comparing each sample to a reference set of calls. This will only have information when providing a validation callset.
* `BATCH-caller.vcf` -- Variants called for a population/batch of samples by a particular caller.
* `BATCH-caller.db` -- A [GEMINI database](https://github.com/arq5x/gemini) associating variant calls with a wide variety of third party annotations. This provides a queryable framework for assessing variant quality statistics.

Sample directories:
* `SAMPLE-caller.vcf` -- Variants calls for an individual sample.
* `SAMPLE-gdc-viral-completeness.txt` -- Optional viral contamination estimates. File is of the format **depth, 1x, 5x, 25x**. **depth** is the number of reads aligning to the virus. **1x, 5x, 25x** are percentage of the viral sequence covered by reads of 1x, 5x, 25x depth. Real viral contamination will have broad coverage across the entire genome, so high numbers for these values, depending on sequencing depth. High depth and low viral sequence coverage means a likely false positive.

### RNA-seq

Project directory:
```
├── annotated_combined.counts -- gene counts with symbols from featureCounts (don't use this)
├── bcbio-nextgen-commands.log -- commands run by bcbio
├── bcbio-nextgen.log -- logging information from bcbio run
├── combined.counts -- gene counts with gene IDs from featureCounts (don't use this)
├── metadata.csv -- provided metadata about each sample
├── multiqc
    ├── multiqc_report.html -- multiQC report
├── programs.txt -- program versions of tools run
├── project-summary.yaml -- YAML description of project, with derived
  metadata
└── tx2gene.csv -- transcript to gene mappings for use with tximport
```
Sample directories:
```
S1
├── S1-ready.bam -- coordinate-sorted whole genome alignments
├── S1-ready.counts -- featureCounts counts (don't use this)
├── S1-transcriptome.bam -- alignments to the transcriptome
├── salmon
│   ├── abundance.h5 -- h5 object, usable with sleuth
│   └── quant.sf -- salmon quantifications, usable with tximport
└── STAR
    ├── S1-SJ.bed -- STAR junction file in BED format
    └── S1-SJ.tab -- STAR junction file in tabular format
```
bcbioRNASeq directory:
```
bcbioRNASeq/
├── data
│ ├── bcb.rda -- bcbioRNASeq object with gene-level data
├── data-transcript
│ ├── bcb.rda -- bcbioRNASeq object with transcript-level data
├── quality_control.html -- quality control report
├── quality_control.Rmd -- RMarkdown that generated quality control report
├── results
│ └── 2019-11-21
│ ├── gene -- gene level information
│ │ └── counts
│ │ ├── counts.csv.gz -- count matrix from tximport, suitable for count-based analyses
│ │ └── metadata.csv.gz -- metadata and quality control data for samples
│ ├── quality_control
│ │ └── tpm.csv.gz -- TPM from tximport, use for visualization
│ └── transcript -- transcript level information
│ └── counts
│ ├── counts.csv.gz -- transcript level count matrix, suitable for count-based analyses needed transcript-level data
│ └── metadata.csv.gz -- metadata and quality control for samples
```
Workflow for analysis:

For gene-level analyses, we recommend loading the gene-level counts.csv.gz and the metadata.csv.gz and using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) to do the analysis. For a more in-depth walkthrough of how to use DESeq2, refer to our [DGE_workshop](https://hbctraining.github.io/DGE_workshop_salmon/schedule/).

For transcript-level analyses, we recommend using [sleuth](https://seqcluster.readthedocs.io/mirna_annotation.html) with the bootstrap samples. You can load the abundance.h5 files from Salmon, or if you set `kallisto` as an expression caller, use the abundance.h5 files from that.

Another great alternative is to use the Salmon quantification to look at differential transcript usage (DTU) instead of differential transcript expression (DTE). The idea behind DTU is you are looking for transcripts of genes that have been flipped from one isoform to another. The [Swimming downstream](https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html#salmon-quantification) tutorial has a nice walkthrough of how to do that.

### single cell RNA-Seq

Project directory:
* `tagcounts.mtx` -- count matrix compatible with dgCMatrix type in R.
* `tagcounts-dupes.mtx` -- count matrix compatible with dgCMatrix type in R but with the duplicated reads counted.
* `tagcounts.mtx.colnames` -- cell names that would be the columns for the matrix.
* `tagcounts.mtx.rownames` -- gene names that would be the rows for the matrix.
* `tagcounts.mtx.metadata` -- metadata that match the colnames for the matrix. This is coming from the barcode.csv file and the metadata given in the YAML config file. for the matrix.
* `cb-histogram.txt` -- total number of dedup reads assigned to a cell. Comparing colSums(tagcounts.mtx) to this number can tell you how many reads mapped to genes.

To create Seurat object:

in bash:
```shell
mkdir data
cd data
result_dir=bcbio_project/final/project_dir
cp $result_dir/tagcounts.mtx matrix.mtx
cp $result_dir/tagcounts.mtx.colnames barcodes.tsv
cp $result_dir/tagcounts.mtx.rownames features.tsv
for f in *;do gzip $f;done;
cd ..
```
in R:
```r
library(Seurat)
counts <- Read10X(data.dir = "data", gene.column = 1)
seurat_object <- CreateSeuratObject(counts = counts, min.features = 100)
saveRDS(seurat_object, "seurat.bcbio.RDS")
```

Sample directories:
* `SAMPLE-transcriptome.bam` -- BAM file aligned to transcriptome.
* `SAMPLE-mtx.*` -- gene counts as explained in the project directory.

### small RNA-seq

Project directory:
* `counts_mirna.tsv` -- miRBase miRNA count matrix.
* `counts.tsv` -- miRBase isomiRs count matrix. The ID is made of 5 tags: miRNA name, SNPs, additions, trimming at 5 and trimming at 3. Here there is detail explanation of the [naming](https://seqcluster.readthedocs.io/mirna_annotation.html).
* `counts_mirna_novel.tsv` -- miRDeep2 miRNA count matrix.
* `counts_novel.tsv` -- miRDeep2 isomiRs. See counts.tsv explanation for more detail. count matrix.
* `seqcluster` -- output of [seqcluster](https://github.com/lpantano/seqcluster) tool. Inside this folder, counts.tsv has count matrix for all clusters found over the genome.
* `seqclusterViz` -- input file for interactive browser at <https://github.com/lpantano/seqclusterViz>
* `report` -- Rmd template to help with downstream analysis like QC metrics, differential expression, and clustering.

Sample directories:
* `SAMPLE-mirbase-ready.counts` -- counts for miRBase miRNAs.
* `SAMPLE-novel-ready` -- counts for miRDeep2 novel miRNAs.
* `tRNA` -- output for [tdrmapper](https://github.com/sararselitsky/tDRmapper).

### ATAC-seq

Sample directories:

Below is an example sample directory for a sample called _rep1_. There are four sets of peak files for each sample, for each peak caller, one set for each of the nucleosome-free (NF), mononucleosome (MF), dinucleosome (DF) and trinucleosome (TF) regions. There are BAM files of reach of those regions as well.
```
rep1
├── macs2 -- peak calls from macs2
│   ├── rep1-NF_control_lambda.bdg.gz -- local lambda estimate for poisson distribution from control samples in bedgraph format, NF regions only
│   ├── rep1-NF_peaks.narrowpeak -- peaks in narrowPeak format in NF regions
│   ├── rep1-NF_summits.bed -- top of peak in bed format of NF regions
│   └── rep1-NF_treat_pileup.bdg.gz -- bedgraph for rep1 sample of NF regions
├── rep1-NF.bam -- BAM of just nucleosome free regions
├── rep1-MN.bam -- BAM of just mononucleosome regions
├── rep1-DN.bam -- BAM of just dinucleosome regions
├── rep1-TN.bam -- BAM of just trinucleosome regions
├── rep1-ready.bam -- BAM of all alignments
└── rep1-ready.bw -- bigwig file of all alignments
```

### Downstream analysis

This section collects useful scripts and tools to do downstream analysis of bcbio-nextgen outputs. If you have pointers to useful tools, please add them to the documentation.

* [Calculate and plot coverage](https://github.com/bcbio/bcbio-nextgen/issues/195#issuecomment-39071048) with matplolib, from Luca Beltrame.
* [Another way](https://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html) to visualize coverage for targeted NGS (exome) experiments with bedtools and R, from Stephen Turner
* assess the efficiency of targeted enrichment sequencing with [ngscat](http://ngscat.clinbioinfosspa.es/start)
