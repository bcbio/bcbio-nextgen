Getting started
---------------

.. _example-pipelines:

Example pipelines
=================

We supply example input configuration files for comparison purposes
and to help in understanding the pipeline.

Exome with validation against reference materials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example calls variants on NA12878 exomes from `EdgeBio's`_
clinical sequencing pipeline, and compares them against
reference materials from NIST's `Genome in a Bottle`_
initiative. This supplies a full regression pipeline to ensure
consistency of calling between releases and updates of third party
software.

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
    $ wget https://s3.amazonaws.com/bcbio_nextgen/NA12878-nist-v2_13-NGv3-pass.vcf.gz
    $ wget https://s3.amazonaws.com/bcbio_nextgen/NA12878-nist-v2_13-NGv3-regions.bed.gz
    $ gunzip NA12878-nist-*.gz
    $ wget https://s3.amazonaws.com/bcbio_nextgen/NGv3.bed.gz
    $ gunzip NGv3.bed.gz

Finally run the analysis, distributed on 8 local cores, with::

    $ mkdir work && cd work
    $ bcbio_nextgen.py bcbio_system.yaml ../input ../config/NA12878-exome-methodcmp.yaml -n 8

The ``grading-summary.csv`` contains detailed comparisons of the results
to the NIST reference materials.

Whole genome
~~~~~~~~~~~~
An input configuration for running whole gnome variant calling with
bwa and GATK, using Illumina's `Platinum genomes project`_
(`NA12878-illumina.yaml`_). To run the analysis:

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
    bcbio_nextgen.py /path/to/your/bcbio_system.yaml ../input ../config/NA12878-illumina.yaml -n 16

- Examine summary of concordance and discordance to comparison calls
  from the ``grading-summary.csv`` file in the work directory.

.. _EdgeBio's: http://www.edgebio.com/
.. _Platinum genomes project: http://www.illumina.com/platinumgenomes/
.. _NA12878-illumina.yaml: https://raw.github.com/chapmanb/bcbio-nextgen/master/config/examples/NA12878-illumina.yaml

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
