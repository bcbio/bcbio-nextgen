Configuration
-------------

Two configuration files, in easy to write `YAML format`_, specify
details about your system and samples to run:

- ``bcbio_system.yaml`` High level information about the system,
  including locations of installed programs like Picard and GATK.
  These apply across multiple runs. The automated installer creates
  a ready to go system configuration file.

- ``bcbio_sample.yaml`` Details about a set of samples to process,
  including input files and analysis options. You configure these for
  each set of samples to process.

Commented `system`_ and `sample`_ example files are available in the
``config`` directory. The :ref:`example-pipelines` section contains
additional examples of ready to run sample files.

Sample information
~~~~~~~~~~~~~~~~~~

The sample configuration file defines ``details`` of each sample to process::

    details:
      - analysis: variant2
        algorithm:
        metadata:
          batch: Batch1
        description: Example1
        genome_build: hg19

- ``analysis`` Analysis method to use [variant2, RNA-seq]
- ``algorithm`` Parameters to configure algorithm inputs. Options
  described in more detail below.
- ``metadata`` Additional descriptive metadata about the sample. The
  ``batch`` input defines a batch that the sample falls in. We perform
  multi-sample variant calling on all samples with the same batch name.
- ``description`` Unique name for this sample. Required.
- ``genome_build`` Genome build to align to, which references a genome
  keyword in Galaxy to find location build files.

Upload
~~~~~~

The ``upload`` section of the sample configuration file defines a
method to extract the final pipeline outputs::

     upload:
       dir: /local/filesystem/directory

General parameters:

- ``method`` Upload method to employ. Defaults to local filesystem.
  [filesystem, galaxy, s3]
- ``dir`` Local filesystem directory to copy to.

Galaxy parameters:

- ``galaxy_url`` URL of the Galaxy instance to upload to. Upload
  assumes you are able to access a shared directory also present on
  the Galaxy machine.
- ``galaxy_api_key`` User API key to access Galaxy: see the
  `Galaxy API`_ documentation.
- ``galaxy_library`` Name of the Galaxy Data Library to upload to. You
  can specify this globally for a project in ``upload`` or for
  individual samples in the sample details section.
- ``galaxy_role`` Specific Galaxy access roles to assign to the
  uploaded datasets. This is optional and will default to the access
  of the parent data library if not supplied. You can specify this
  globally for a project in ``upload`` or for individual samples in
  the sample details section.

S3 parameters:

- ``bucket`` AWS bucket to upload to
- ``access_key_id`` AWS access key ID from Amazon credentials page
- ``secret_access_key`` AWS secret key ID from Amazon credentials page
- ``reduced_redundancy`` Flag to determine if we should store S3 data
  with reduced redundancy: cheaper but less reliable [false, true]

Algorithm parameters
~~~~~~~~~~~~~~~~~~~~

The YAML configuration file provides a number of hooks to customize
analysis in the sample configuration file. Place these under the
``algorithm`` keyword.

Alignment
=========

-  ``aligner`` Aligner to use: [bwa, bowtie, bowtie2, mosaik, novoalign,
   false]
-  ``bam_clean`` Clean an input BAM when skipping alignment step. This
   handles adding read groups, sorting to a reference genome and
   filtering problem records that cause problems with GATK. Set to
   ``picard`` to do Picard/GATK based cleaning.
-  ``bam_sort`` Allow sorting of input BAMs when skipping alignment
   step (``aligner`` set to false). Options are coordinate or
   queryname. For additional processing through standard pipelines
   requires coordinate sorted inputs. The default is to not do
   additional sorting and assume pre-sorted BAMs.
-  ``trim_reads`` Whether to trim off 3' B-only ends from fastq reads
   [false, true]
-  ``align_split_size``: Split FASTQ files into specified number of
   records per file. Allows parallelization at the cost of increased
   temporary disk space usage.
-  ``quality_bin``: Perform binning of quality scores with CRAM to
   reduce file sizes. Uses the Illumina 8-bin approach. Supply a list
   of times to perform binning: [prealignment, postrecal]
-  ``quality_format`` Quality format of fastq inputs [illumina,
   standard]
-  ``write_summary`` Write a PDF summary of results [true, false]
-  ``merge_bamprep`` Merge regional BAM prepped files into a final
   prepared BAM. false avoids the time consuming merge when you only
   want variant calls [true, false]
-  ``coverage_bigwig`` Generate a bigwig file of coverage, for loading
   into the UCSC genome browser [true, false]

Experimental information
========================

-  ``coverage_interval`` Regions covered by sequencing. Influences GATK
   options for filtering [exome, genome, regional]
-  ``coverage_depth`` Depth of sequencing coverage. Influences GATK
   variant calling [high, low]
-  ``hybrid_target`` BED file with target regions for hybrid selection
   experiments.
-  ``ploidy`` Ploidy of called reads. Defaults to 2 (diploid).

Variant calling
===============

-  ``variantcaller`` Variant calling algorithm. Can be a list of
   multiple options [gatk, freebayes, varscan, samtools,
   gatk-haplotype, cortex]
-  ``validate`` A VCF file of expected variant calls to perform
    validation and grading of output variants from the pipeline.
    This provides a mechanism to ensure consistency of calls against
    a known set of variants, supporting comparisons to genotyping
    array data or reference materials.
- ``validate_regions`` A BED file of regions to evaluate in. This
  defines specific regions covered by the ``validate`` VCF  file.
- ``validate_genome_build``: Genome build of the validation file, if
  different than the samples genome build. Helps manage hg19/GRCh37
  chromosome naming differences.
-  ``variant_regions`` BED file of regions to call variants in.
-  ``mark_duplicates`` Identify and remove variants [false, true]
-  ``recalibrate`` Perform variant recalibration [true, false]
-  ``realign`` Type of variant alignment to perform, Defaults to
   GATK realignment. [gatk, gkno, false]
-  ``phasing`` Do post-call haplotype phasing of variants. Defaults to
   no phasing [false, gatk]

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

Parallelization
===============

- ``nomap_split_size`` Unmapped base pair regions required to split
  analysis into blocks. Creates islands of mapped reads surrounded by
  unmapped (or N) regions, allowing each mapped region to run in
  parallel. (default: 100)

- ``nomap_split_targets`` Number of target intervals to attempt to
  split processing into. This picks unmapped regions evenly spaced
  across the genome to process concurrently. Limiting targets prevents
  a large number of small targets. (default: 2000)

Ensemble variant calling
========================

In addition to single method variant calling, we support calling with
multiple calling methods and consolidating into a final Ensemble
callset. This requires the `bcbio.variation`_ toolkit to perform the
consolidation. An example configuration in the ``algorithm`` section is::

    variantcaller: [gatk, freebayes, samtools, gatk-haplotype, varscan]
    ensemble:
      format-filters: [DP < 4]
      classifier-params:
        type: svm
      classifiers: 
        balance: [AD, FS, Entropy]
        calling: [ReadPosEndDist, PL, PLratio, Entropy, NBQ]
      trusted-pct: 0.65

The ``ensemble`` set of parameters configure how to combine calls from
the multiple methods:

- ``format-filters`` A set of filters to apply to variants before
  combining. The example removes all calls with a depth of less than
  4.
- ``classifier-params`` Parameters to configure the machine learning
  approaches used to consolidate calls. The example defines an SVM
  classifier.
- ``classifiers`` Groups of classifiers to use for training and
  evaluating during machine learning. The example defines two set of
  criteria for distinguishing reads with allele balance issues and
  those with low calling support.
- ``trusted-pct`` Define threshold of variants to include in final
  callset. In the example, variants called by more than 65% of the
  approaches (4 or more callers) pass without being requiring SVM
  filtering.

.. _config-resources:
   
Resources
~~~~~~~~~

The ``resources`` section allows customization of locations of programs
and memory and compute resources to devote to them::

    resources:
      bwa:
        cores: 12
        cmd: /an/alternative/path/to/bwa
      gatk:
        jvm_opts: ["-Xms2g", "-Xmx4g"]
        dir: /usr/share/java/gatk

- ``cmd`` Location of an executable. By default, we assume executables
  are on the path.
- ``dir`` For software not distributed as a single executable, like
  files of Java jars, the location of the base directory.
- ``cores`` Cores to use for multi-proccessor enabled software. On
  cluster systems, match this with the number of physical cores
  available on individual machines.
- ``jvm_opts`` Specific memory usage options for Java software.

Resources will continue to expand to allow direct customization of
commandline options as well as fine grained control over research
usage.

.. _bcbio.variation: https://github.com/chapmanb/bcbio.variation
.. _CloudBioLinux: https://github.com/chapmanb/cloudbiolinux
.. _YAML format: https://en.wikipedia.org/wiki/YAML#Examples
.. _GATK resource bundle: http://www.broadinstitute.org/gsa/wiki/index.php/GATK_resource_bundle
.. _GATK: http://www.broadinstitute.org/gatk/
.. _system: https://github.com/chapmanb/bcbio-nextgen/blob/master/config/bcbio_system.yaml
.. _sample: https://github.com/chapmanb/bcbio-nextgen/blob/master/config/bcbio_sample.yaml
.. _Galaxy API: http://wiki.galaxyproject.org/Learn/API
