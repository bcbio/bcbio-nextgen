Configuration
-------------

Two configuration files, in easy to write `YAML format`_, specify
details about your system and samples to run:

- ``bcbio_system.yaml`` High level information about the system,
  including locations of installed programs like Picard and GATK.
  These apply across multiple runs.

- ``bcbio_sample.yaml`` Details about a set of samples to process,
  including input files and analysis options. You configure these for
  each set of samples to process.

Commented example files are available in the ``config`` directory:

- `example system config`_
- `example sample config`_

Options
~~~~~~~

The YAML configuration file provides a number of hooks to customize
analysis in the sample configuration file. Place these under the
``analysis`` keyword. For variant calling:

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

Broad's `GATK`_ pipeline drives variant (SNP and Indel) analysis.
This requires some associated data files, and also has some configurable
options. The relevant section from the ``bcbio_system.yaml`` file is::

    dbsnp: variation/dbsnp_132.vcf
    train_hapmap: variation/hapmap_3.3.vcf
    train_1000g_omni: variation/1000G_omni2.5.vcf
    train_indels: variation/Mills_Devine_2hit.indels.vcf

The dbSNP and training files are from the `GATK resource bundle`_. These
are inputs into the training models for recalibration. The automated
`CloudBioLinux`_ data scripts will download and install these in the
variation subdirectory relative to the genome files.

.. _CloudBioLinux: https://github.com/chapmanb/cloudbiolinux
.. _YAML format: https://en.wikipedia.org/wiki/YAML#Examples
.. _GATK resource bundle: http://www.broadinstitute.org/gsa/wiki/index.php/GATK_resource_bundle
.. _GATK: http://www.broadinstitute.org/gatk/
.. _example system config: https://github.com/chapmanb/bcbb/blob/master/nextgen/config/bcbio_system.yaml
.. _example sample config: https://github.com/chapmanb/bcbb/blob/master/nextgen/config/bcbio_sample.yaml

