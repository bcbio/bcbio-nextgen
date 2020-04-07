## Getting started

This example calls variants using NA12878 exome data from [EdgeBio's](https://www.edgebio.com/) clinical sequencing pipeline, and compares them against reference materials from NIST's [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) initiative.

1. Install bcbio python package and tools
    ```shell
    wget https://raw.github.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
    python bcbio_nextgen_install.py [bcbio_installation_path] --tooldir=[tools_installation_path] --nodata
    ```

2. Install hg38 reference genome and bwa indices
    ```shell
    bcbio_nextgen.py update -u skip --genomes hg38 --aligners bwa
    ```
    See more detailed instructions in the installation user story.

3. Get the input configuration file, fastq reads, reference materials and analysis regions:
    ```shell
    mkdir -p NA12878-exome-eval
    cd NA12878-exome-eval
    wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-exome-methodcmp-getdata.sh
    bash NA12878-exome-methodcmp-getdata.sh
    ```

4.  Run the analysis, distributed on 8 local cores, with:
    ```shell
    cd work
    bcbio_nextgen.py ../config/NA12878-exome-methodcmp.yaml -n 8
    ```
    The `grading-summary.csv` contains detailed comparisons of the results to the NIST reference materials, enabling rapid comparisons of methods.

### Project structure

bcbio encourages a project structure:
```
project/
├── config
├── final
├── input
└── work
```
with the `project.yaml` configuration in the `config` directory, the input files (fastq, bam, bed) in the `input` directory, the outputs of the pipeline in the `final` directory, and the actual processing done in the `work` directory.

Typical bcbio run:
* copy or link input files in the `input` directory
* set pipeline parameters in `config/project.yaml`
* run the `bcbio_nextgen.py` script from inside the `work`
* review the results in `final`
* delete `work` with intermediate files.

### Logging

There are 3 logging files in the `log` directory within your working folder:
* `bcbio-nextgen.log` High level logging information about the analysis. This provides an overview of major processing steps and useful checkpoints for assessing run times.
* `bcbio-nextgen-debug.log` Detailed information about processes including stdout/stderr from third party software and error traces for failures. Look here to identify the status of running pipelines or to debug errors. It labels each line with the hostname of the machine it ran on to ease debugging in distributed cluster environments.
* `bcbio-nextgen-commands.log` Full command lines for all third party software tools run.

### Example pipelines

We supply example input configuration files for validation and to help in understanding the pipeline.

#### Whole genome trio (50x) - hg38

This input configuration runs whole genome bwa alignment and GATK variant calling. It uses a father/mother/child trio from the [CEPH NA12878 family](https://blog.goldenhelix.com/wp-content/uploads/2013/03/Utah-Pedigree-1463-with-NA12878.png): NA12891, NA12892, NA12878. Illumina's [Platinum genomes project](https://www.illumina.com/platinumgenomes.html) has 50X whole genome sequencing of the three members. The analysis compares results against a reference NA12878 callset from NIST's [Genome in a
Bottle](https://www.nist.gov/programs-projects/genome-bottle) initiative.

To run the analysis do:
```shell
mkdir NA12878-trio-eval
cd NA12878-trio-eval
mkdir config input work
cd config
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-trio-wgs-validate.yaml
cd ../input
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-trio-wgs-validate-getdata.sh
bash NA12878-trio-wgs-validate-getdata.sh
cd ../work
bcbio_nextgen.py ../config/NA12878-trio-wgs-validate.yaml -n 16
```
This is a large whole genome analysis and meant to test both pipeline scaling and validation across the entire genome. It can take multiple days to run depending on available cores. It requires 300GB for the input files and 1.3TB for the work directory. Smaller examples below exercise the pipeline with less disk and computational requirements.

#### Cancer tumor normal - GRCh37

This example calls variants using multiple approaches in a paired tumor/normal cancer sample from the [ICGC-TCGA DREAM challenge](https://www.synapse.org/#!Synapse:syn312572/wiki/58893). It uses [synthetic dataset 3](https://www.synapse.org/#!Synapse:syn312572/wiki/62018) which has multiple subclones, enabling detection of lower frequency variants. Since the dataset is freely available and has a truth set, this allows us to do a full evaluation of variant callers.

To get the data:
```shell
mkdir -p cancer-dream-syn3/config cancer-dream-syn3/input cancer-dream-syn3/work
cd cancer-dream-syn3/config
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/cancer-dream-syn3.yaml
cd ../input
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/cancer-dream-syn3-getdata.sh
bash cancer-dream-syn3-getdata.sh
```
Run with:
```shell
cd ../work
bcbio_nextgen.py ../config/cancer-dream-syn3.yaml -n 8
```
The configuration and data file has downloads for exome only and whole genome analyses. It enables exome by default, but you can use the larger whole genome evaluation by uncommenting the relevant parts of the configuration and retrieval script.

#### Cancer-like mixture with Genome in a Bottle samples

This example simulates somatic cancer calling using a mixture of two Genome in a Bottle samples, NA12878 as the "tumor" mixed with NA24385 as the background. The [Hartwig Medical Foundation](https://www.hartwigmedicalfoundation.nl/en/) and [Utrecht Medical Center](https://www.umcutrecht.nl/en/Research/Strategic-themes/Cancer) generated this "tumor/normal" pair by physical mixing of samples prior to sequencing. The GiaB FTP directory has [more details on the design and truth sets](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/use_cases/mixtures/UMCUTRECHT_NA12878_NA24385_mixture_10052016/README-NA12878_NA24385_mixture.txt). The sample has variants at 15% and 30%, providing the ability to look at lower frequency mutations.

To get the data:
```shell
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/cancer-giab-na12878-na24385-getdata.sh
bash cancer-giab-na12878-na24385-getdata.sh
```
Then run the analysis with:
```shell
cd work
bcbio_nextgen.py ../config/cancer-giab-na12878-na24385.yaml -n 16
```

#### Structural variant calling -- whole genome NA12878 (50x)

This example runs structural variant calling with multiple callers (Lumpy, Manta and CNVkit), providing a combined output summary file and validation metrics against NA12878 deletions. It uses the same NA12878 input as the whole genome trio example.

To run the analysis do:
```shell
mkdir -p NA12878-sv-eval
cd NA12878-sv-eval
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-sv-getdata.sh
bash NA12878-sv-getdata.sh
cd work
bcbio_nextgen.py ../config/NA12878-sv.yaml -n 16
```
This is large whole genome analysis and the timing and disk space requirements for the NA12878 trio analysis above apply here as well.

#### Human genome build 38

Validate variant calling on human genome build 38, using two different builds (with and without alternative alleles) and three different validation datasets (Genome in a Bottle prepared with two methods and Illumina platinum genomes). To run:
```shell
mkdir -p NA12878-hg38-val
cd NA12878-hg38-val
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-hg38-validate-getdata.sh
bash NA12878-hg38-validate-getdata.sh
mkdir -p work
cd work
bcbio_nextgen.py ../config/NA12878-hg38-validate.yaml -n 16
```

#### Whole genome (10x)

An input configuration for running whole gnome variant calling with bwa
and GATK, using Illumina's [Platinum genomes project](https://www.illumina.com/platinumgenomes.html) ([NA12878-illumina.yaml](https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-illumina.yaml)). See this [blog post on whole genome scaling](https://bcb.io/2013/05/22/scaling-variant-detection-pipelines-for-whole-genome-sequencing-analysis/) for expected run times and more information about the pipeline. To run the analysis:

* Create an input directory structure like:
    ```shell
    ├── config
    │   └── NA12878-illumina.yaml
    ├── input
    └── work
    ```
* Retrieve inputs and comparison calls:
    ```shell
    cd input
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR091/ERR091571/ERR091571_1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR091/ERR091571/ERR091571_2.fastq.gz
    ```
* Retrieve configuration input file:
    ```shell
    cd config
    wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-illumina.yaml
    ```
* Run analysis on 16 core machine:
    ```shell
    cd work
    bcbio_nextgen.py ../config/NA12878-illumina.yaml -n 16
    ```
* Examine summary of concordance and discordance to comparison calls from the `grading-summary.csv` file in the work directory.
