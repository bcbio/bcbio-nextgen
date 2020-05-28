# Getting started

## Workflow

This example calls variants using NA12878 exome data from
[EdgeBio's](https://www.edgebio.com/) clinical sequencing pipeline,
and compares them against reference materials from NIST's
[Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) initiative.

### 1. Install bcbio python package and tools
    ```shell
    wget https://raw.github.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
    python bcbio_nextgen_install.py [bcbio_installation_path] --tooldir=[tools_installation_path] --nodata
    ```

### 2. Install hg38 reference genome and bwa indices
    ```shell
    bcbio_nextgen.py update -u skip --genomes hg38 --aligners bwa
    ```
    See more detailed instructions in the installation user story.

### 3. Get the input configuration file, fastq reads, reference materials and analysis regions:
    ```shell
    mkdir -p NA12878-exome-eval
    cd NA12878-exome-eval
    wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-exome-methodcmp-getdata.sh
    bash NA12878-exome-methodcmp-getdata.sh
    ```

### 4. Run the analysis, distributed on 8 local cores, with:
    ```shell
    cd work
    bcbio_nextgen.py ../config/NA12878-exome-methodcmp.yaml -n 8
    ```
    Parameters of the analysis are specified in the yaml configuration file:
    ```
    upload:
      dir: ../final
    details:
      - files: [../input/NA12878-NGv3-LAB1360-A_1.fastq.gz, ../input/NA12878-NGv3-LAB1360-A_2.fastq.gz]
      description: NA12878
      metadata:
        sex: female
      analysis: variant2
      genome_build: hg38
      algorithm:
        aligner: bwa
        variantcaller: gatk-haplotype
        validate: giab-NA12878/truth_small_variants.vcf.gz
        validate_regions: giab-NA12878/truth_regions.bed
        variant_regions: capture_regions/Exome-NGv3
    ```
    Running time is ~2h.

### 5. Explore results in `NA12878-exome-eval/final`:
*  `date_project/multiqc` - quality contol
*  `date_project/NA12878-gatk-haplotype-annotated.vcf.gz` - annotated variants
*  `NA12878/NA12878-callable.bed` - callable regions
*  `final/NA12878/NA12878-ready.bam` - bam file
*  `date_project/bcbio-nextgen-commands.log` - commands ran to produce results
*  `date_project/grading-summary-NA12878.csv` - validation results. 
False Discovery Rate (FDR) for SNPs here is 3% (i.e. 97% precision for SNPs), 
so the precision is quite low. 
One reason of low precision could be that NA12878-NGv3-LAB1360 WES dataset
was sequenced in 2013 or earlier, so it could be of somewhat lower quality.
We left it here for educational purpose. 
With a modern NA12878 dataset you can achieve >99% precision and >99% sensitivity using bcbio/gatk, 
see [germline variants user story](germline_variants.html#workflow1-validate-hg38-calls).
Comparing QC and validations in the two NA12878 WES datasets illustrates how sequencing quality affects variant calling precision and sensitivity. 
Another point one could make when comparing the two validations is that NA12878-NGv3-LAB1360 
has a larger target (133,288 SNPs vs 37,033), so the choice of `variant_regions` directly influences validation results.
Including only regions with high coverage, excluding low complexity regions leads to increased precision.
A larger bed file with more regions included is a more stressful test for combination of capture kit/sequencing instrument/aligner/variant caller/filters.

## What is next?
Bcbio documentation is organized by user stories. We support 22 user stories (extended use cases):
* 14 data processing user stories corresponding to different types of NGS data
and biological questions
* 8 infrastructural stories.

### Data processing stories
1. Somatic variants
2. Bulk RNA-seq expression
3. Single cell RNA-seq
4. HLA typing
5. Germline small variants
6. 3'prime digital gene expression
7. Structural variants
8. ChIP/ATAC-seq
9. Methylation
10. Bulk RNA-seq variants
11. Bulk RNA-seq fusion
12. Fast RNA-seq
13. Disambiguation
14. Small RNA-seq

### Infrastructural stories
1. Getting started
2. Installation
3. Configuration
4. Parallel execution
5. Outputs
6. Development
7. Cloud
8. CWL

### A typical user story contains:
- `workflow` - step-by step description of how to run the pipeline with example data
- `parameters` - describes yaml config file parameters relevant for the user story
- `output` - describes output files
- `steps` - outlines low level step the pipeline performs to process data
- `validation` - validation results when available
- `description`
- `references`

Try running bcbio with your own data, report [issues](https://github.com/bcbio/bcbio-nextgen/issues),
contribute to the codebase and documentation on [github](https://github.com/bcbio/bcbio-nextgen/).
