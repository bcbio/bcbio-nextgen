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
batched population level variant calls.

Project directory
~~~~~~~~~~~~~~~~~
- ``project-summary.yaml`` -- Top level YAML format summary file with
  statistics on read alignments and duplications as well as analysis
  specific metrics.
- ``programs.txt`` -- Program versions for bcbio-nextgen and software
  run in the pipeline. This enables reproduction of analyses.

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

Downstream analysis
===================

This section collects useful scripts and tools to do downstream analysis of
bcbio-nextgen outputs. If you have pointers to useful tools, please add them to
the documentation.

- `Calculate and plot coverage`_ with matplolib, from Luca Beltrame.

.. _Calculate and plot coverage:  https://github.com/chapmanb/bcbio-nextgen/issues/195#issuecomment-39071048
