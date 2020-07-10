# Bulk RNA-seq

## Workflow

This example aligns and creates count files for use with downstream analyses using a subset of the SEQC data from the [FDA's Sequencing Quality Control project](http://www.fdaseqc.org/).

### 1. Install STAR index
```
bcbio_nextgen.py upgrade -u skip --genomes hg38 --aligners star --cores 10
```

### 2. Setup bcbio project
Download input data, create project structure and config files. This will download six samples from the SEQC project, three from the HBRR panel and three from the UHRR panel (100G download). :
```shell
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/rnaseq-seqc-getdata.sh
bash rnaseq-seqc-getdata.sh
```

Step 2 in detail:

#### 2.1 Create input directory and download fastq files
```
mkdir -p input
```

#### 2.2 Download a template yaml file describing RNA-seq analysis
```
wget --no-check-certificate https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/rnaseq-seqc.yaml
```
rnaseq-seqc.yaml:
```
# Template for human RNA-seq using Illumina prepared samples
---
details:
  - analysis: RNA-seq
    genome_build: hg38
    algorithm:
      aligner: star
      strandedness: unstranded
upload:
  dir: ../final
```

#### 2.3 Download fastq files into input dir
```
cd input
for SAMPLE in SRR950078 SRR950079 SRR950080 SRR950081 SRR950082 SRR950083
do 
   wget -c -O ${SAMPLE}_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/${SAMPLE}/${SAMPLE}_1.fastq.gz
   wget -c -O ${SAMPLE}_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/${SAMPLE}/${SAMPLE}_2.fastq.gz
done
```

#### 2.4 Prepare a sample sheet
```
wget -c --no-check-certificate https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/seqc.csv
```

seqc.csv:
```
samplename,description,panel
SRR950078,UHRR_rep1,UHRR
SRR950079,HBRR_rep1,HBRR
SRR950080,UHRR_rep2,UHRR
SRR950081,HBRR_rep2,HBRR
SRR950082,UHRR_rep3,UHRR
SRR950083,HBRR_rep3,HBRR
SRR950084,UHRR_rep4,UHRR
```

#### 2.5 Generate yaml config file for analysis
project.yaml = template.yaml x sample_sheet.csv
```
cd ../
bcbio_nextgen.py -w template rnaseq-seqc.yaml input/seqc.csv input
```

In the result you should see a folder structure:
```
seqc
|---config
|---final
|---work
```

`seqc/config/seqc.yaml` is the main config file to run the bcbio project.

### 3. Run the analysis:
```shell
cd seqc/work
bcbio_nextgen.py ../config/seqc.yaml -n 8
```

This will run a full scale RNAseq experiment using STAR as the aligner and will take a long time to finish on a single machine. At the end it will output counts, Cufflinks quantitation and a set of QC results about each lane. If you have a cluster you can [parallelize it](parallel) to speed it up considerably.

## Parameters

* `transcript_assembler` If set, will assemble novel genes and transcripts and merge the results into the known annotation. Can have multiple values set in a list. Supports ['cufflinks', 'stringtie'].
* `transcriptome_align` If set to True, will also align reads to just the transcriptome, for use with EBSeq and others.
* `expression_caller` A list of optional expression callers to turn on. Supports ['cufflinks', 'express', 'stringtie', 'dexseq', 'kallisto']. Salmon and count based expression estimation are run by default. Sailfish is deprecated.
* `fusion_caller` A list of optional fusion callers to turn on. Supports [oncofuse, pizzly, ericscript, arriba].
* `variantcaller` Variant calling algorithm to call variants on RNA-seq data. Supports [gatk-haplotype] or [vardict].
* `spikein_fasta` A FASTA file of spike in sequences to quantitate.
* `quantify_genome_alignments` If set to True, run Salmon quantification using the genome alignments from STAR, when available. If STAR alignments are not available, use Salmon's SA mode with decoys.
* `bcbiornaseq` A dictionary of key-value pairs to be passed as options to bcbioRNAseq. Currently supports _organism_ as a key and takes the latin name of the genome used (_mus
musculus_, _homo sapiens_, etc) and _interesting_groups_ which will be used to color
quality control plots:
    ```yaml
    bcbiornaseq:
      organism: homo sapiens
      interesting_groups: [treatment, genotype, etc, etc]
    ```
You will need to also turn on `bcbiornaseq` by turning it on via `tools_on: [bcbiornaseq]`.

## Output

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

## Description

RNA-seq pipeline includes steps for quality control, adapter trimming, alignment, variant calling, transcriptome reconstruction and post-alignment quantitation at the level of the gene and isoform.

We recommend using the STAR aligner for all genomes. Use Tophat2 only if you do not have enough RAM available to run STAR (about 30 GB).

Our current recommendation is to run adapter trimming only if using the Tophat2 aligner. Adapter trimming is very slow, and aligners that soft clip the ends of reads such as STAR and hisat2, or algorithms using pseudoalignments like Salmon handle contaminant sequences at the ends properly. This makes trimming unnecessary. Tophat2 does not perform soft clipping so if using Tophat2, trimming must still be done.

Salmon, which is an extremely fast alignment-free method of quantitation, is run for all experiments. Salmon can accurately quantitate the expression of genes, even ones which are hard to quantitate with other methods (see [this paper](https://doi.org/10.1186/s13059-015-0734-x) for example for Sailfish, which performs similarly to Salmon). Salmon can also quantitate at the transcript level which can help gene-level analyses
(see [this paper](https://doi.org/10.12688/f1000research.7563.1) for example). We recommend using the Salmon quantitation rather than the counts from featureCounts to perform downstream quantification.

Although we do not recommend using the featureCounts based counts, the alignments are still useful because they give you many more quality metrics than the quasi-alignments from Salmon.

After a bcbio RNA-seq run there will be in the `upload` directory a directory for each sample which contains a BAM file of the aligned and unaligned reads, a `salmon` directory with the output of Salmon, including TPM values, and a `qc` directory with plots from FastQC and qualimap.

In addition to directories for each sample, in the `upload` directory there is a project directory which contains a YAML file describing some summary statistics for each sample and some provenance data about the bcbio run. In that directory is also a `combined.counts` file with the featureCounts derived counts per cell.
