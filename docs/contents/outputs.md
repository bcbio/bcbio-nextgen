# Outputs

bcbio-nextgen runs in a temporary work directory which contains a number of processing intermediates. Pipeline completion extracts the final useful output files into a separate directory, specified by the [Upload](contents/configuration:upload). This configuration allows upload to local directories, Galaxy, or Amazon S3. Once extracting and confirming the output files, you can delete the temporary directory to save space.

The output directory contains sample specific output files labeled by sample name and a more general project directory. The sample directories contain all of the sample specific output files, while the project directory contains global files like project summaries or batched population level variant calls. See the [Teaching](teaching) documentation for a full variant calling example with additional details about configuration setting and resulting output files.

## Project directory:
* `project-summary.yaml` -- Top level YAML format summary file with statistics on read alignments and duplications as well as analysis specific metrics.
* `programs.txt` -- Program versions for bcbio-nextgen and software run in the pipeline. This enables reproduction of analyses.
* `multiqc` -- [MultiQC](https://multiqc.info/) report. multiqc_report.html combines quality metrics from multiple tools (listed in multiqc_config.yaml).
  - General statistics/Dup - % of duplicates among mapped reads, calculated by [bcbio.qc.qualimap.py](https://github.com/bcbio/bcbio-nextgen/blob/69bc24d703d3a0166caabf833ddd9e514ff1d445/bcbio/qc/qualimap.py#L219)
  - General statistics/% Dups - duplication statistics from [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). It uses first 200k reads to generate % DUP, not a good proxy for RNA-seq data, [read more in fastqc docs](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html)
  - Fastqc/Sequence counts/duplicates - same as General statistics/%Dups
  - `General statistics/ontarget_pct` = 100.0 * ontarget / mapped_unique, see [code](https://github.com/bcbio/bcbio-nextgen/blob/a3473775db06540c10b5f20ddc2043b8cc99d1f8/bcbio/qc/coverage.py#L63)
  - `General statistics/Usable_pct` = 100.0 * ontarget / total_reads
* `metadata.csv` -- CSV with the metadata in the YAML file.
* `data_versions.csv` -- Data versions for bcbio-nextgen and software

## Sample directories:
* `SAMPLE/qc` -- Directory of quality control runs for the sample. These include charts and metrics for assessing quality of sequencing and analysis.
* `SAMPLE-ready.bam` -- A prepared BAM file of the aligned reads. Depending on the analysis used, this may include trimmed, recalibrated and realigned reads following alignment.

## Why do I have so many coverage metrics? Which one should I use?

## Interpretation of ontarget_pct vs usable_pct
Usually, `ontarget_pct` > `usable_pct` for WES without UMI because `mapped_unique < total_reads`, `ontarget_pct ~ usable_pct` for projects with UMI, because UMI consensus reads are mapped by definition.
Sometimes `ontarget_pct` < `usable_pct`, the difference is <=0.2% for UMI projects, because `total_reads` is calculated by samtools, and `mapped_unique` by [readstats.py](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/bam/readstats.py#L37), which leads to a slight difference in `mapped_unique > total`.

## Interpretation of mosdepth median coverage vs qualimap median coverage
`Mosdepth median coverage` < `Qualimap median coverage`. For example, mosdepth = 133, qualimap = 169. Mosdepth calculates the median over the average coverages of contigs (chromosomes), see `mosdepth/Average coverage per contig` figure. Most of the chromosomes have coverage 150-200X, but chromosomes are very few and additional chromosomes have coverage 20-30X, so the median is lowered by that fact. qualimap goes for a median over all exonic regions, see `Qualimap/Coverage Histogram` figure. Exonic regions are many and the median is higher. Use qualimap median coverage as a better estimate. Note that `tools_on: qualimap_full` should be used to get reliable estimates from qualimap. If you are subsetting bam for qualimap (default), use mosdepth median coverage.

## Interpretation of bcbio(mosdepth) average target coverage vs qualimap mean coverage
First of all, these two should not be confused with median coverage from mosdepth and qualimap (mean vs median).
`bcbio_average_target` is calculated by [get_average_coverage](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/variation/coverage.py#L154), in the end it is coverage from `mosdepth`: `NA12878-exome-eval/work/coverage/NA12878/NA12878-variant_regions.regions.bed.gz ` for NA12878 WES project. It is calculated **excluding** duplicated reads and reads with unmapped mate. See also statistics here in the `bcbio project/work/coverage/mapped_stats.txt`. Qualimap includes these reads, so usually for WES data with and without UMIs `bcbio_average_target < qualimap_mean_coverage`. However, we've seen some projects with UMIs and panels where `bcbio_average_target_coverage >> qualimap_mean_coverage`. Use `bcbio_average`.

## Downstream analysis

This section collects useful scripts and tools to do downstream analysis of bcbio-nextgen outputs. If you have pointers to useful tools, please add them to the documentation.

* [Calculate and plot coverage](https://github.com/bcbio/bcbio-nextgen/issues/195#issuecomment-39071048) with matplolib, from Luca Beltrame.
* [Another way](https://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html) to visualize coverage for targeted NGS (exome) experiments with bedtools and R, from Stephen Turner
* assess the efficiency of targeted enrichment sequencing with [ngscat](http://ngscat.clinbioinfosspa.es/start)
