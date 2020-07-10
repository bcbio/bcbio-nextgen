# smallRNA-seq

## Overview
bcbio supports configurable best-practices pipeline for smallRNA-seq quality controls,
adapter trimming, miRNA/isomiR quantification and other small RNA detection.

* Adapter trimming:
  * [atropos](https://atropos.readthedocs.io/en/latest/guide.html)
  * [dnapi](https://github.com/jnktsj/DNApi) for adapter de-novo detection
* Sequence alignment:
  * [STAR](https://code.google.com/archive/p/rna-star) for genome annotation
  * bowtie, _bowtie2_ and [hisat2](https://daehwankimlab.github.io/hisat2/) for genome annotation as an option
* Specific small RNAs quantification (miRNA/tRNAs...):
  * [seqbuster](https://github.com/lpantano/seqbuster) for miRNA annotation
  * [MINTmap](https://github.com/TJU-CMC-Org/MINTmap) for tRNA fragments annotation
  * [miRge2](https://github.com/mhalushka/miRge) for alternative small RNA quantification. To setup this tool, you need to install manually miRge2.0, and download the library data for your species. Read how to install and download the [data](https://github.com/mhalushka/miRge#download-libraries). If you have `human` folder at `/mnt/data/human` the option to pass to resources will be `/mnt/data`. Then setup `resources`:
    ```yaml
    resources:
        mirge:
            options: ["-lib $PATH_TO_PARENT_SPECIES_LIB"]
    ```
* Quality control: [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* Other small RNAs quantification:
  * [seqcluster](https://github.com/lpantano/seqcluster)
  * [mirDeep2](https://www.mdc-berlin.de/content/mirdeep2-documentation) for miRNA prediction

The pipeline generates a _RMD_ template file inside `report` folder that can be rendered with knitr. An example of the report is [here](https://github.com/lpantano/mypubs/blob/master/srnaseq/mirqc/ready_report.md). Count table (`counts_mirna.tst`) from mirbase miRNAs will be inside `mirbase` or final project folder. Input files for [isomiRs](https://github.com/lpantano/isomiRs) package for isomiRs analysis will be inside each sample in `mirbase` folder. If mirdeep2 can run, count table (`counts_mirna_novel.tsv`) for novel miRNAs will be inside `mirdeep2` or final project folder. tdrmapper results will be inside each sample inside `tdrmapper` or final project folder.

## Parameters
* `adapters` The 3' end adapter that needs to be remove. For NextFlex protocol you can add `adapters: ["4N", "$3PRIME_ADAPTER"]`. For any other options you can use resources: `atropos:options:["-u 4", "-u -4"]`.
* `species` 3 letters code to indicate the species in mirbase classification (i.e. hsa for human).
* `aligner` Currently STAR is the only one tested although bowtie can be used as well.
* `expression_caller` A list of expression callers to turn on: trna, seqcluster, mirdeep2, mirge
* `transcriptome_gtf` An optional GTF file of the transcriptome to for seqcluster.
* `spikein_fasta` A FASTA file of spike in sequences to quantitate.
* `umi_type: 'qiagen_smallRNA_umi'` Support of Qiagen UMI small RNAseq protocol.

## Output
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
