# ATAC-seq
The ATAC-seq pipeline in bcbio follows recommendations from the [ENCODE ATAC-seq pipeline](https://www.encodeproject.org/atac-seq/) and [Yiwei Niu's excellent guide](https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/). 

bcbio aligns the reads and cleans up the alignments, removing duplicates,
multimappers and reads aligning to mitochondria. It then breaks up the BAM files
into separate BAM files for nucleosome free (NF), mononucleosome (MN),
dinucleosome (DN) and trinucleosome (TN) fractions and calls peaks separately on
each fraction and also calls peaks on all of the fractions together.

Consensus peaks of the nucleosome free peaks are created by choosing the peak
with the highest score when peaks overlap, described more in depth in the
[bedops
documentation](https://bedops.readthedocs.io/en/latest/content/usage-examples/master-list.html).
A matrix of peak counts is created with featureCounts that can be used with
downstream count-based differential expression callers like DESeq2/limma/edgeR.

ENCODE quality control metrics and other quality control information is added to the [MultiQC](https://multiqc.info) report. A separate ATAC-seq specific quality control report is generated using [ataqv](https://github.com/ParkerLab/ataqv).

## Description of example dataset
We will be using [ENCSR312LQX](https://www.encodeproject.org/experiments/ENCSR312LQX) and
[ENCSR310MLB](https://www.encodeproject.org/experiments/ENCSR310MLB/) from the ENCODE project
as our example datasets. These are P0 samples from the mouse hindbrain and mouse forebrain, each
with a single replicate each. We'll use these samples to call differential affinity between the mouse 
hindbrain and forebrain.

### 1. Download the example data and configuration files
This downloads the input data, creates the project structure and example configuration files.

#### 1.1 Create input directory and download FASTQ files.
```bash
mkdir atac-example
cd atac-example
mkdir -p fastq

# hindbrain samples
wget --no-check-certificate https://www.encodeproject.org/files/ENCFF547YID/@@download/ENCFF547YID.fastq.gz -O fastq/hindbrain_rep1_R1.fastq.gz
wget --no-check-certificate https://www.encodeproject.org/files/ENCFF131VHT/@@download/ENCFF131VHT.fastq.gz -O fastq/hindbrain_rep1_R2.fastq.gz
wget --no-check-certificate https://www.encodeproject.org/files/ENCFF971XEA/@@download/ENCFF971XEA.fastq.gz -O fastq/hindbrain_rep2_R1.fastq.gz
wget --no-check-certificate https://www.encodeproject.org/files/ENCFF215WAD/@@download/ENCFF215WAD.fastq.gz -O fastq/hindbrain_rep2_R2.fastq.gz
# forebrain samples
wget --no-check-certificate https://www.encodeproject.org/files/ENCFF296GZG/@@download/ENCFF296GZG.fastq.gz -O fastq/forebrain_rep1_R1.fastq.gz
wget --no-check-certificate https://www.encodeproject.org/files/ENCFF664RZO/@@download/ENCFF664RZO.fastq.gz -O fastq/forebrain_rep1_R2.fastq.gz
wget --no-check-certificate https://www.encodeproject.org/files/ENCFF197GTC/@@download/ENCFF197GTC.fastq.gz -O fastq/forebrain_rep2_R1.fastq.gz
wget --no-check-certificate https://www.encodeproject.org/files/ENCFF209GGJ/@@download/ENCFF209GGJ.fastq.gz -O fastq/forebrain_rep2_R2.fastq.gz
```

#### 1.2 Download template YAML file describing the ATAC-seq analysis

```bash
mkdir -p metadata
wget --no-check-certificate http://s3.amazonaws.com/bcbio-nextgen/atac_userstory_data/atac-example.yaml -O metadata/atac-example.yaml
```

atac-example.yaml: 
```yaml
details:
- analysis: chip-seq
  genome_build: mm10
  algorithm:
    aligner: bwa
    peakcaller: [macs2]
    chip_method: atac
    keep_duplicates: False
    keep_multimapped: False
upload:
  dir: ../final
```

#### 1.3 Create a sample sheet

```bash
wget --no-check-certificate http://s3.amazonaws.com/bcbio-nextgen/atac_userstory_data/hindbrain_forebrain.csv -O metadata/hindbrain_forebrain.csv
```

hindbrain_forebrain.csv:
```
samplename,description,region,replicate
forebrain_rep1_R1.fastq.gz,forebrain_rep1,forebrain,rep1
forebrain_rep2_R1.fastq.gz,forebrain_rep2,forebrain,rep2
hindbrain_rep1_R1.fastq.gz,hindbrain_rep1,hindbrain,rep1
hindbrain_rep2_R1.fastq.gz,hindbrain_rep2,hindbrain,rep2
```

The only two fields required in this file are `samplename` and `description`, you can put whatever you want
for the other columns. We recommend adding any additional metdata you know about the samples here.

### 2. Generate YAML config file for analysis
```bash
bcbio_nextgen.py -w template metadata/atac-example.yaml metadata/hindbrain_forebrain.csv fastq
```

In the result you should see a folder structure:
```
hindbrain_forebrain
|---config
|---final
|---work
```

`hindbrain_forebrain/config/hindbrain_forebrain.yaml` is the main config file to run the bcbio project. You will
see this file has a copy of the parameters in `atac-example.yaml` for each sample.

### 3. Run the analysis
This will run the analysis on a local machine, using 16 cores.
```bash
cd hindbrain_forebrain/work
bcbio_nextgen.py ../config/hindbrain_forebrain.yaml -n 16
```

## Parameters
* `peakcaller`: `[macs2]` bcbio just supports MACS2
* `aligner`: supports `bowtie2` and `bwa`. `bwa` will result in a superset of the peaks called by `bowtie2`.
* `chip_method`: set to `atac` to run the ATAC-seq pipeline
* `keep_duplicates`: do not remove duplicates before peak calling. Defaults to _False_.
* `keep_multimapped`: do not remove multimappers before peak calling. Defaults to _False_.

## Output

### Project directory

```
├── 2020-05-01_hindbrain_forebrain
│   ├── ataqv
│   │   ├── index.html -- QC report from ataqv
│   ├── bcbio-nextgen-commands.log -- list of commands run by bcbio
│   ├── bcbio-nextgen.log -- stdout of bcbio log
│   ├── consensus 
│   │   ├── consensus.bed -- consensus peaks from NF fraction
│   │   ├── consensus-counts.tsv -- table of alignments per peak for each sample, calculated by featureCounts
│   ├── data_versions.csv -- versions of data used in the pipeline
│   ├── metadata.csv -- supplied metadata about the samples
│   ├── multiqc
│   │   ├── multiqc_report.html -- multiQC report with useful quality control metrics
│   ├── programs.txt -- versions of programs run in the pipeline
```

### Sample directories

```
├── forebrain_rep1
│   ├── forebrain_rep1-DN.bam -- dinucleosome alignments
│   ├── forebrain_rep1-full.bam -- all fraction alignments
│   ├── forebrain_rep1-MN.bam -- mononucleosome alignments
│   ├── forebrain_rep1-NF.bam -- nucleosome-free alignments
│   ├── forebrain_rep1-ready.bam -- identifical to -full
│   ├── forebrain_rep1-ready.bam.bai 
│   ├── forebrain_rep1-ready.bw -- bigwig file of full alignments
│   ├── forebrain_rep1-TN.bam -- trinucleosome alignments
│   ├── macs2 -- contains peak calls for each fraction, including the full peak calls
```

## Downstream analysis

### Quality Control
The **MultiQC** report in the project directory under `multiqc/multiqc_report.html`
and at the **ataqv** report in the project directory under
`ataqv/ataqv_report.html` have useful quality control information that you can
use to help decide if your ATAC-seq project worked.

It is hard to give specific cutoffs of metrics to use since the kit, the sample
material, the organism, the genome annotations and so on all affect all of the
metrics. We generally look at the samples as a whole for an experiment and see
if any of the samples are outliers in the important metrics. In the **MultiQC**
report, we look at the percentage of reads in the peaks, the mapping percentage,
the 
[ENCODE library complexity statistics](https://www.encodeproject.org/data-standards/terms/) and the FastQC
metrics to try to spot samples with problems.

In the **ataqv** report, we look at the HQAA fragment length distribution plot.
Ideally, this plot should show a periodic uptick every 200 bases, which
corresponds to the different nucleosome fractions. The samples should be
enriched for < 100 which is the nucleosome free fraction, 200 for the
mononucleosome fraction, 400 for the dinucleosome fraction and 600 for the
trinucleosome fraction. Often you will not see this behavior though even in
libraries that were successful. But if some of your samples have this and others
do not, that is something to be concerned about.

You should see an enrichment around the transcription start sites, if you are
missing that then your experiment likely failed. The **peaks** table in the
**tables** tab in the **ataqv** report has a measurement of the high quality
autosomal alignments overlapping peaks, **ataqv** calculates this metric using
all of the peaks, not just the peaks from the nucleosome-free fraction, so this
is useful to look at as well. See the [ataqv github
repository](https://github.com/ParkerLab/ataqv/issues/13) for a discussion of
the ranges of values you can expect to see for metrics in the **ataqv** report
along with other values to look at that might be informative. The Parker lab
reprocessed samples from many publications with **ataqv** and posted the reports
[here](https://theparkerlab.med.umich.edu/data/porchard/ataqv-public-survey/)
which is helpful to browse through to get an idea of what ranges of values you
can expect. As you can see, they can be all over the place.

#### hindbrain vs forebrain QC reports
- [MultiQC report](http://atac-userstory.s3-website.us-east-2.amazonaws.com/multiqc_report.html)
- [ataqv report](http://atac-userstory.s3-website.us-east-2.amazonaws.com)

### Differential affinity analysis
For doing differential affinity analysis we recommend loading the consensus peak
table from the `consensus/consensus.counts` file in the project directory. This
is a table of counts per peak in the nucleosome-free fraction for each sample
that you can use in any standard count-based differential expression tools like
[DESeq2](
https://bioconductor.org/packages/release/bioc/html/DESeq2.html)/[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)/[limma](https://bioconductor.org/packages/release/bioc/html/limma.html).
Often you will find tutorials using
[DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html)
but that uses these callers under the hood, so you can just call them directly
and skip an intermediate step if you want. Either way works. The DiffBind tutorials
are great for understanding how to go about with your downstream analyses.

#### hindbrain vs forebrain differential affinity reports
- [RMarkdown](http://atac-userstory.s3-website.us-east-2.amazonaws.com/peaks.Rmd)
- [HTML report](http://atac-userstory.s3-website.us-east-2.amazonaws.com/peaks.html)
- [example data](http://atac-userstory.s3-website.us-east-2.amazonaws.com/differential-affinity-example.tar.gz)
