Configuration
-------------

Analysis
~~~~~~~~~~~~~~~~~~~~

The YAML configuration file provides a number of hooks to customize
analysis. Place these under the ``analysis`` keyword. For variant
calling:

-  ``aligner`` Aligner to use: [bwa, bowtie, bowtie2, mosaik, novoalign,
   false]
-  ``trim_reads`` Whether to trim off 3' B-only ends from fastq reads
   [false, true]
-  ``align_split_size``: Split FASTQ files into specified number of
   records per file. Allows parallelization at the cost of increased
   temporary disk space usage.
-  ``variantcaller`` Variant calling algorithm [gatk, freebayes]
-  ``quality_format`` Quality format of fastq inputs [illumina,
   standard]
-  ``coverage_interval`` Regions covered by sequencing. Influences GATK
   options for filtering [exome, genome, regional]
-  ``coverage_depth`` Depth of sequencing coverage. Influences GATK
   variant calling [high, low]
-  ``hybrid_target`` BED file with target regions for hybrid selection
   experiments.
-  ``variant_regions`` BED file of regions to call variants in.
-  ``ploidy`` Ploidy of called reads. Defaults to 2 (diploid).
-  ``recalibrate`` Perform variant recalibration [true, false]
-  ``mark_duplicates`` Identify and remove variants [false, true]
-  ``realign`` Do variant realignment [true, false]
-  ``write_summary`` Write a PDF summary of results [true, false]

Global reference files for variant calling and assessment:

-  ``train_hapmap``, ``train_1000g_omni``, ``train_indels`` Training
   files for GATK variant recalibration.
-  ``call_background`` Background VCF to use for calling.

Variant calling
~~~~~~~~~~~~~~~

Broad's `GATK`_ pipeline drives variant (SNPs and Indels) analysis.
This requires some associated data files, and also has some configurable
options. The relevant section from the ``post_process.yaml`` file is:

::

    coverage_depth: "low" # other options: high
    coverage_interval: "exome" # other options: genome, regional
    dbsnp: variation/dbsnp_132.vcf
    train_hapmap: variation/hapmap_3.3.vcf
    train_1000g_omni: variation/1000G_omni2.5.vcf
    train_indels: variation/Mills_Devine_2hit.indels.vcf

The dbSNP and training files are from the `GATK resource bundle`_. These
are inputs into the training models for recalibration. The automated
[CloudBioLinux][o6] data scripts will download and install these in the
variation subdirectory relative to the genome files.

The ``coverage_depth`` and ``coverage_interval`` are adjustable from the
defaults within individual ``run_info.yaml`` files describing the
sample; a fastq file with standard phred quality scores, full genome
coverage and high sequencing depth:

::

    - analysis: SNP calling
      algorithm:
        quality_format: Standard
        coverage_interval: genome
        coverage_depth: high

.. _GATK resource bundle: http://www.broadinstitute.org/gsa/wiki/index.php/GATK_resource_bundle
.. _GATK: http://www.broadinstitute.org/gatk/

