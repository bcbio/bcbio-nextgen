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
    wget https://s3.amazonaws.com/bcbio_nextgen/NA12878-illumina-example.vcf.gz
    gunzip NA12878-illumina-example.vcf.gz

- Retrieve configuration input file::

    cd config
    wget https://raw.github.com/chapmanb/bcbio-nextgen/master/config/examples/NA12878-illumina.yaml

- Run analysis on 16 core machine::

    cd work
    bcbio_nextgen.py ../input ../config/NA12878-illumina.yaml -n 16

- Examine summary of concordance and discordance to comparison calls
  from the ``grading-summary.csv`` file in the work directory.

.. _EdgeBio's: http://www.edgebio.com/
.. _Platinum genomes project: http://www.illumina.com/platinumgenomes/
.. _NA12878-illumina.yaml: https://raw.github.com/chapmanb/bcbio-nextgen/master/config/examples/NA12878-illumina.yaml
.. _blog post on whole genome scaling: http://bcbio.wordpress.com/2013/05/22/scaling-variant-detection-pipelines-for-whole-genome-sequencing-analysis/

Exome with validation against reference materials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example calls variants on NA12878 exomes from `EdgeBio's`_
clinical sequencing pipeline, and compares them against
reference materials from NIST's `Genome in a Bottle`_
initiative. This supplies a full regression pipeline to ensure
consistency of calling between releases and updates of third party
software. This is a large full exome examples with multiple aligners
and variant callers, so can take more than 24 hours on machines using
multiple cores.

First get the input configuration file::

    $ mkdir config && cd config
    $ wget https://raw.github.com/chapmanb/bcbio-nextgen/master/config/\
       examples/NA12878-exome-methodcmp.yaml

Then the fastq reads, reference materials and analysis regions::

    $ cd .. && mkdir input && cd input
    $ wget https://dm.genomespace.org/datamanager/file/Home/EdgeBio/\
       CLIA_Examples/NA12878-NGv3-LAB1360-A/NA12878-NGv3-LAB1360-A_1.fastq.gz
    $ wget https://dm.genomespace.org/datamanager/file/Home/EdgeBio/\
       CLIA_Examples/NA12878-NGv3-LAB1360-A/NA12878-NGv3-LAB1360-A_2.fastq.gz
    $ wget https://s3.amazonaws.com/bcbio_nextgen/NGv3.bed.gz
    $ wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/variant_calls/NIST/NIST_IntegratedCalls_12datasets_130517_HetHomVarPASS_VQSRv2.15.vcf.gz
    $ wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/variant_calls/NIST/union12callableMQonlymerged_addcert_nouncert_excludesegdups_excludedecoy_noCNVs_v2.15b.bed.gz
    $ gunzip *.vcf.gz
    $ gunzip *.bed.gz

Finally run the analysis, distributed on 8 local cores, with::

    $ mkdir work && cd work
    $ bcbio_nextgen.py ../config/NA12878-exome-methodcmp.yaml -n 8

The ``grading-summary.csv`` contains detailed comparisons of the results
to the NIST reference materials.

Note that this example requires a full licensed version of novoalign,
since it uses gzipped inputs and multicore processing. You can still
run the example with only bwa alignments by removing lanes 1 and 3
from the ``NA12878-exome-methodcmp.yaml`` sample configuration file.

Exome with Ensemble calling
~~~~~~~~~~~~~~~~~~~~~~~~~~~

An example configuration for running `ensemble variant calling`_ on
multiple exome samples (`NA12878-ensemble.yaml`_).

.. _NA12878-ensemble.yaml: https://raw.github.com/chapmanb/bcbio-nextgen/master/config/examples/NA12878-ensemble.yaml
.. _ensemble variant calling: http://bcbio.wordpress.com/2013/02/06/an-automated-ensemble-method-for-combining-and-evaluating-genomic-variants-from-multiple-callers/
.. _Genome in a Bottle: http://www.genomeinabottle.org/

Test suite
==========

The test suite exercises the scripts driving the analysis, so are a good
starting point to ensure correct installation. Run tests from the main
code directory using `nose`_. To test the main variant calling
pipeline::

     $ cd tests
     $ nosetests -v -s -a speed=1

To run the full test suite::

     $ nosetest -v -s

``tests/test_automated_analysis.py`` exercises the full framework using
an automatically downloaded test dataset. It runs through barcode
deconvolution, alignment and full SNP analysis. Tweak the configuration
for the tests for your environment:

-  ``tests/data/automated/post_process.yaml`` -- May need adjustment to
   point to installed software in non-standard locations. Change the
   num\_cores parameter to test multiple processor and parallel
   execution.
-  ``tests/data/automated/run_info.yaml`` -- Change the ``analysis``
   variable can to 'Standard' if variant calling is not required in your
   environment. This will run a smaller pipeline of alignment and
   analysis.

.. _nose: http://somethingaboutorange.com/mrl/projects/nose/
