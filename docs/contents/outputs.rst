Outputs
-------
bcbio-nextgen runs in a temporary work directory which contains a
number of processing intermediates. Pipeline completion extracts the
final useful output files into a separate directory, specified by the
:ref:`upload-configuration`. This configuration allows upload to local
directories, Galaxy, or Amazon S3. Once extracting and confirming the
output files, you can delete the temporary directory to save space.

Common files
============

The output directory contains sample specific output files labeled by
sample name and a more general project directory. The sample
directories contain all of the sample specific output files, while the
project directory contains global files like project summaries or
batched population level variant calls. See the :ref:`teaching` documentation
for a full variant calling example with additional details about configuration
setting and resulting output files.

Project directory
~~~~~~~~~~~~~~~~~
- ``project-summary.yaml`` -- Top level YAML format summary file with
  statistics on read alignments and duplications as well as analysis
  specific metrics.
- ``programs.txt`` -- Program versions for bcbio-nextgen and software
  run in the pipeline. This enables reproduction of analyses.
- ``multiqc`` run `MultiQC`_ to gather all QC metrics from different tools, such as,
  cutadapt, featureCounts, samtools, STAR ... into an unique HTML report.
- ``metadata.csv`` -- CSV with the metadata in the YAML file.

.. _MultiQC: http://multiqc.info

Sample directories
~~~~~~~~~~~~~~~~~~
- ``SAMPLE/qc`` -- Directory of quality control runs for the sample.
  These include charts and metrics for assessing quality of sequencing
  and analysis.
- ``SAMPLE-ready.bam`` -- A prepared BAM file of the aligned reads.
  Depending on the analysis used, this may include trimmed,
  recalibrated and realigned reads following alignment.

Variant calling
===============

Project directory
~~~~~~~~~~~~~~~~~

- ``grading-summary.csv`` -- Grading details comparing each sample to
  a reference set of calls. This will only have information when
  providing a validation callset.
- ``BATCH-caller.vcf`` -- Variants called for a population/batch of
  samples by a particular caller.
- ``BATCH-caller.db`` -- A `GEMINI database`_ associating variant
  calls with a wide variety of third party annotations. This provides
  a queryable framework for assessing variant quality statistics.

.. _GEMINI database: https://github.com/arq5x/gemini

Sample directories
~~~~~~~~~~~~~~~~~~
- ``SAMPLE-caller.vcf`` -- Variants calls for an individual sample.
- ``SAMPLE-gdc-viral-completeness.txt`` -- Optional viral contamination estimates. File is of the format **depth, 1x, 5x, 25x**. **depth** is the number of reads aligning to the virus. **1x, 5x, 25x** are percentage of the viral sequence covered by reads of 1x, 5x, 25x depth. Real viral contamination will have broad coverage across the entire genome, so high numbers for these values, depending on sequencing depth. High depth and low viral sequence coverage means a likely false positive.

RNA-seq
=======

Project directory
~~~~~~~~~~~~~~~~~

- ``annotated_combined.counts`` -- featureCounts counts matrix
  with gene symbol as an extra column.
- ``combined.counts`` -- featureCounts counts matrix
  with gene symbol as an extra column.
- ``combined.dexseq`` -- DEXseq counts matrix with
  exonID as first column.
- ``combined.gene.sf.tmp`` -- Sailfish gene count
  matrix normalized to TPM.
- ``combined.isoform.sf.tpm`` -- Sailfish transcript
  count matix normalized to TPM.
- ``combined.sf`` -- Sailfish raw output, all samples
  files are pasted one after another.
- ``tx2gene.csv`` -- Annotation file needed for DESeq2
  to use Sailfish output.

Sample directories
~~~~~~~~~~~~~~~~~~

- ``SAMPLE-transcriptome.bam`` -- BAM file aligned to transcriptome.
- ``SAMPLE-ready.counts`` -- featureCounts gene counts output.
- ``salmon`` -- Salmon output.

single cell RNA-seq
===================

Project directory
~~~~~~~~~~~~~~~~~

- ``tagcounts.mtx`` -- count matrix compatible with dgCMatrix type in R.
- ``tagcounts-dupes.mtx`` -- count matrix compatible with dgCMatrix type in R
  but with the duplicated reads counted.
- ``tagcounts.mtx.colnames`` -- cell names that would be the columns
  for the matrix.
- ``tagcounts.mtx.rownames`` -- gene names that would be the rows
  for the matrix.
- ``tagcounts.mtx.metadata`` -- metadata that match the colnames
  for the matrix. This is coming from the barcode.csv file and
  the metadata given in the YAML config file.
  for the matrix.
- ``cb-histogram.txt`` -- total number of dedup reads assigned to a cell.
  Comparing colSums(tagcounts.mtx) to this number can tell you how many
  reads mapped to genes.

To create Seurat object:

in bash::

  mkdir data
  cd data
  result_dir=bcbio_project/final/project_dir
  cp $result_dir/tagcounts.mtx matrix.mtx
  cp $result_dir/tagcounts.mtx.colnames barcodes.tsv
  cp $result_dir/tagcounts.mtx.rownames features.tsv
  for f in *;do gzip $f;done;
  cd ..

in R::

  library(Seurat)
  counts <- Read10X(data.dir = "data", gene.column = 1)
  seurat_object <- CreateSeuratObject(counts = counts, min.features = 100)
  saveRDS(seurat_object, "seurat.bcbio.RDS")

Sample directories
~~~~~~~~~~~~~~~~~~

- ``SAMPLE-transcriptome.bam`` -- BAM file aligned to transcriptome.
- ``SAMPLE-mtx.*`` -- gene counts as explained in the project directory.


small RNA-seq
=============

Project directory
~~~~~~~~~~~~~~~~~

- ``counts_mirna.tsv`` -- miRBase miRNA
  count matrix.
- ``counts.tsv`` -- miRBase isomiRs count matrix. The ID is made of 5 tags:
  miRNA name, SNPs, additions, trimming at 5 and trimming at 3.
  Here there is detail explanation of the `naming`_ .
- ``counts_mirna_novel.tsv`` -- miRDeep2 miRNA
  count matrix.
- ``counts_novel.tsv`` -- miRDeep2 isomiRs. See counts.tsv explanation for more detail.
  count matrix.
- ``seqcluster`` -- output of `seqcluster`_ tool.
  Inside this folder, counts.tsv has count matrix
  for all clusters found over the genome.
- ``seqclusterViz`` -- input file for interactive
  browser at https://github.com/lpantano/seqclusterViz
- ``report`` -- Rmd template to help with downstream
  analysis like QC metrics, differential expression, and
  clustering.

.. _naming: http://seqcluster.readthedocs.io/mirna_annotation.html

Sample directories
~~~~~~~~~~~~~~~~~~

- ``SAMPLE-mirbase-ready.counts`` -- counts for miRBase miRNAs.
- ``SAMPLE-novel-ready`` -- counts for miRDeep2 novel miRNAs.
- ``tRNA`` -- output for `tdrmapper`_.

.. _seqcluster: https://github.com/lpantano/seqcluster
.. _tdrmapper: https://github.com/sararselitsky/tDRmapper

Downstream analysis
===================

This section collects useful scripts and tools to do downstream analysis of
bcbio-nextgen outputs. If you have pointers to useful tools, please add them to
the documentation.

- `Calculate and plot coverage`_ with matplolib, from Luca Beltrame.
- `Another way`_ to visualize coverage for targeted NGS (exome) experiments with bedtools and R, from Stephen Turner
- assess the efficiency of targeted enrichment sequencing with `ngscat`_

.. _ngscat: http://www.bioinfomgp.org/ngscat
.. _Calculate and plot coverage:  https://github.com/bcbio/bcbio-nextgen/issues/195#issuecomment-39071048
.. _Another way: http://gettinggeneticsdone.blogspot.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html
