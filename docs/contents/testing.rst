Getting started
---------------

Overview
========

1. Create a `sample configuration file`_ for your project
   (substitute the example BAM and fastq names below with the full
   path to your sample files)::

         bcbio_nextgen.py -w template gatk-variant project1 sample1.bam sample2_1.fq sample2_2.fq

   This uses a standard template (GATK best practice variant calling)
   to automate creation of a full configuration for all samples. See
   :ref:`automated-sample-config` for more details on running the
   script, and manually edit the base template or final output
   file to incorporate project specific configuration. The example
   pipelines provide a good starting point and the
   :ref:`sample-configuration` documentation has full details on
   available options.

2. Run analysis, distributed across 8 local cores::

         bcbio_nextgen.py bcbio_sample.yaml -n 8

3. Read the :ref:`docs-config` documentation for full details on
   adjusting both the sample and system configuration files to match
   your experiment and computational setup.

.. _sample configuration file: https://github.com/chapmanb/bcbio-nextgen/blob/master/config/bcbio_sample.yaml

Logging
=======

There are 3 logging files in the ``log`` directory within your working folder:

- ``bcbio-nextgen.log`` High level logging information about the analysis.
  This provides an overview of major processing steps and useful
  checkpoints for assessing run times.
- ``bcbio-nextgen-debug.log`` Detailed information about processes
  including stdout/stderr from third party software and error traces
  for failures. Look here to identify the status of running pipelines
  or to debug errors. It labels each line with the hostname of the
  machine it ran on to ease debugging in distributed cluster
  environments.
- ``bcbio-nextgen-commands.log`` Full command lines for all third
  party software tools run.

.. _example-pipelines:

Example pipelines
=================

We supply example input configuration files for comparison purposes
and to help in understanding the pipeline.

Whole genome
~~~~~~~~~~~~
An input configuration for running whole gnome variant calling with
bwa and GATK, using Illumina's `Platinum genomes project`_
(`NA12878-illumina.yaml`_). See this
`blog post on whole genome scaling`_ for expected run times and more
information about the pipeline. To run the analysis:

- Create an input directory structure like::

    ├── config
    │   └── NA12878-illumina.yaml
    ├── input
    └── work

- Retrieve inputs and comparison calls::

    cd input
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR091/ERR091571/ERR091571_1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR091/ERR091571/ERR091571_2.fastq.gz
    wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/variant_calls/NIST/\
     NISTIntegratedCalls_13datasets_130719_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.17_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz
    wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/variant_calls/NIST/\
     union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.17.bed.gz
    gunzip *.vcf.gz *.bed.gz

- Retrieve configuration input file::

    cd config
    wget https://raw.github.com/chapmanb/bcbio-nextgen/master/config/examples/NA12878-illumina.yaml

- Run analysis on 16 core machine::

    cd work
    bcbio_nextgen.py ../config/NA12878-illumina.yaml -n 16

- Examine summary of concordance and discordance to comparison calls
  from the ``grading-summary.csv`` file in the work directory.

.. _Platinum genomes project: http://www.illumina.com/platinumgenomes/
.. _NA12878-illumina.yaml: https://raw.github.com/chapmanb/bcbio-nextgen/master/config/examples/NA12878-illumina.yaml
.. _blog post on whole genome scaling: http://bcbio.wordpress.com/2013/05/22/scaling-variant-detection-pipelines-for-whole-genome-sequencing-analysis/

Exome with validation against reference materials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example calls variants on NA12878 exomes from `EdgeBio's`_
clinical sequencing pipeline, and compares them against reference
materials from NIST's `Genome in a Bottle`_ initiative. This supplies
a full regression pipeline to ensure consistency of calling between
releases and updates of third party software. The pipeline performs
alignment with bwa mem and variant calling with FreeBayes, GATK
UnifiedGenotyper and GATK HaplotypeCaller. Finally it integrates all 3
variant calling approaches into a `combined ensemble callset`_.

This is a large full exome example with multiple variant callers, so
can take more than 24 hours on machines using multiple cores.

First get the input configuration file::

    mkdir config && cd config
    wget https://raw.github.com/chapmanb/bcbio-nextgen/master/config/\
     examples/NA12878-exome-methodcmp.yaml

Then the fastq reads, reference materials and analysis regions::

    cd .. && mkdir input && cd input
    wget https://dm.genomespace.org/datamanager/file/Home/EdgeBio/\
     CLIA_Examples/NA12878-NGv3-LAB1360-A/NA12878-NGv3-LAB1360-A_1.fastq.gz
    wget https://dm.genomespace.org/datamanager/file/Home/EdgeBio/\
     CLIA_Examples/NA12878-NGv3-LAB1360-A/NA12878-NGv3-LAB1360-A_2.fastq.gz
    wget https://s3.amazonaws.com/bcbio_nextgen/NGv3.bed.gz
    wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/variant_calls/NIST/\
     NISTIntegratedCalls_13datasets_130719_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.17_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz
    wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/variant_calls/NIST/\
     union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.17.bed.gz
    gunzip *.vcf.gz *.bed.gz

Finally run the analysis, distributed on 8 local cores, with::

    cd .. & mkdir work && cd work
    bcbio_nextgen.py ../config/NA12878-exome-methodcmp.yaml -n 8

The ``grading-summary.csv`` contains detailed comparisons of the results
to the NIST reference materials, enabling rapid comparisons of methods.

.. _combined ensemble callset: http://bcbio.wordpress.com/2013/02/06/an-automated-ensemble-method-for-combining-and-evaluating-genomic-variants-from-multiple-callers/
.. _Genome in a Bottle: http://www.genomeinabottle.org/
.. _EdgeBio's: http://www.edgebio.com/

Cancer tumor normal
~~~~~~~~~~~~~~~~~~~

This example calls variants in a paired cancer sample with tumor/normal
sequencing data. using raw data from `Han et al in PLoS One
<http://www.plosone.org/article/info:doi/10.1371/journal.pone.0064271>`_. This
is a work in progress and we welcome contributions. The goal is to use a full
evaluation dataset to compare calling methods:

Get the input configuration file::

    mkdir config && cd config
    wget https://raw.github.com/chapmanb/bcbio-nextgen/master/config/\
     examples/cancer-paired.yaml

Get fastq reads and analysis regions::

    cd .. && mkdir input && cd input
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR256/ERR256785/ERR256785_1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR256/ERR256785/ERR256785_2.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR256/ERR256786/ERR256786_1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR256/ERR256786/ERR256786_2.fastq.gz
    wget https://gist.github.com/chapmanb/8322238/raw/131a5710ac17039e8e2d350e00a88898e030a958/ERP002442-targeted.bed

Run::

    cd .. & mkdir work && cd work
    bcbio_nextgen.py ../config/cancer-paired.yaml -n 8

Test suite
==========

The test suite exercises the scripts driving the analysis, so are a
good starting point to ensure correct installation. Tests use the
`nose`_ test runner pre-installed as part of the pipeline. Grab the latest
source code::

     $ git clone https://github.com/chapmanb/bcbio-nextgen.git

To run the standard tests::

     $ cd bcbio-nextgen/tests
     $ ./run_tests.sh

To run specific subsets of the tests::

     $ ./run_tests.sh rnaseq
     $ ./run_tests.sh speed=2
     $ ./run_tests.sh devel

By default the test suite will use your installed system configuration
for running tests, substituting the test genome information instead of
using full genomes. If you need a specific testing environment, copy
``tests/data/automated/post_process-sample.yaml`` to
``tests/data/automated/post_process.yaml`` to provide a test-only
configuration.

.. _nose: http://somethingaboutorange.com/mrl/projects/nose/
