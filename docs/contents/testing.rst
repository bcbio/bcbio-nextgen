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

.. _sample configuration file: https://github.com/bcbio/bcbio-nextgen/blob/master/config/bcbio_sample.yaml

Project directory
=================

bcbio encourages a project structure like::

    my-project/
    ├── config
    ├── final
    └── work

with the input configuration in the ``config`` directory, the outputs of the
pipeline in the ``final`` directory, and the actual processing done in the
``work`` directory. Run the ``bcbio_nextgen.py`` script from inside the ``work``
directory to keep all intermediates there.  The ``final`` directory, relative to
the parent directory of the ``work`` directory, is the default location
specified in the example configuration files and gets created during
processing. The ``final`` directory has all of the finished outputs and you can
remove the ``work`` intermediates to cleanup disk space after confirming the
results. All of these locations are configurable and this project structure is
only a recommendation.

.. _logging-output:

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

We supply example input configuration files for validation
and to help in understanding the pipeline.

Whole genome trio (50x)
~~~~~~~~~~~~~~~~~~~~~~~

This input configuration runs whole genome variant calling using bwa, GATK
HaplotypeCaller and FreeBayes. It uses a father/mother/child
trio from the `CEPH NA12878 family`_: NA12891, NA12892, NA12878.
Illumina's `Platinum genomes project`_ has 50X whole genome sequencing of the
three members. The analysis compares results against a reference
NA12878 callset from NIST's `Genome in a Bottle`_ initiative.

To run the analysis do::

  mkdir -p NA12878-trio-eval/config NA12878-trio-eval/input NA12878-trio-eval/work
  cd NA12878-trio-eval/config
  wget https://raw.github.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-trio-wgs-validate.yaml
  cd ../input
  wget https://raw.github.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-trio-wgs-validate-getdata.sh
  bash NA12878-trio-wgs-validate-getdata.sh
  cd ../work
  bcbio_nextgen.py ../config/NA12878-trio-wgs-validate.yaml -n 16

This is a large whole genome analysis and meant to test both pipeline scaling
and validation across the entire genome. It can take multiple days to run
depending on available cores. It requires 300Gb for the input files and 1.3Tb
for the work directory. Smaller examples below exercise the pipeline with
less disk and computational requirements.

.. _CEPH NA12878 family: http://blog.goldenhelix.com/wp-content/uploads/2013/03/Utah-Pedigree-1463-with-NA12878.png

We also have a more extensive evaluation that includes 2 additional variant
callers, Platypus and samtools, and 3 different methods of calling variants:
single sample, pooled, and incremental joint calling. This uses the same input
data as above but a different input configuration file::

  mkdir -p NA12878-trio-eval/work_joint
  cd NA12878-trio-eval/config
  wget https://raw.github.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-trio-wgs-joint.yaml
  cd ../work_joint
  bcbio_nextgen.py ../config/NA12878-trio-wgs-joint.yaml -n 16

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

First get the input configuration file, fastq reads, reference materials and analysis regions::

    mkdir -p NA12878-exome-eval
    cd NA12878-exome-eval
    wget https://raw.github.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-exome-methodcmp-getdata.sh
    bash NA12878-exome-methodcmp-getdata.sh

Finally run the analysis, distributed on 8 local cores, with::

    cd work
    bcbio_nextgen.py ../config/NA12878-exome-methodcmp.yaml -n 8

The ``grading-summary.csv`` contains detailed comparisons of the results
to the NIST reference materials, enabling rapid comparisons of methods.

.. _combined ensemble callset: http://bcb.io/2013/02/06/an-automated-ensemble-method-for-combining-and-evaluating-genomic-variants-from-multiple-callers/
.. _Genome in a Bottle: http://www.genomeinabottle.org/
.. _EdgeBio's: http://www.edgebio.com/

.. _example-cancer:

Cancer tumor normal
~~~~~~~~~~~~~~~~~~~

This example calls variants using multiple approaches in a paired tumor/normal
cancer sample from the `ICGC-TCGA DREAM challenge
<https://www.synapse.org/#!Synapse:syn312572>`_. It uses `synthetic dataset 3
<https://www.synapse.org/#!Synapse:syn312572/wiki/62018>`_ which has multiple
subclones, enabling detection of lower frequency variants. Since the dataset is
freely available and has a truth set, this allows us to do a full evaluation of
variant callers.

To get the data::

    mkdir -p cancer-dream-syn3/config cancer-dream-syn3/input cancer-dream-syn3/work
    cd cancer-dream-syn3/config
    wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/cancer-dream-syn3.yaml
    cd ../input
    wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/cancer-dream-syn3-getdata.sh
    bash cancer-dream-syn3-getdata.sh

Run with::

    cd ../work
    bcbio_nextgen.py ../config/cancer-dream-syn3.yaml -n 8

The configuration and data file has downloads for exome only and whole genome
analyses. It enables exome by default, but you can use the larger whole genome
evaluation by uncommenting the relevant parts of the configuration and retrieval
script.

Cancer-like mixture with Genome in a Bottle samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example simulates somatic cancer calling using a mixture of two Genome in a
Bottle samples, NA12878 as the "tumor" mixed with NA24385 as the background.
The `Hartwig Medical Foundation <http://www.hartwigmedicalfoundation.nl/en/>`_
and `Utrecht Medical Center
<http://www.umcutrecht.nl/en/Research/Research-programs/Cancer>`_ generated this
"tumor/normal" pair by physical mixing of samples prior to sequencing. The GiaB
FTP directory has `more details on the design and truth sets
<ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/use_cases/mixtures/UMCUTRECHT_NA12878_NA24385_mixture_10052016/README-NA12878_NA24385_mixture.txt>`_.
The sample has variants at 15% and 30%, providing the ability to look at lower
frequency mutations.

To get the data::

    wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/cancer-giab-na12878-na24385-getdata.sh
    bash cancer-giab-na12878-na24385-getdata.sh

Then run the analysis with::

    cd work
    bcbio_nextgen.py ../config/cancer-giab-na12878-na24385.yaml -n 16

Structural variant calling -- whole genome NA12878 (50x)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example runs structural variant calling with multiple callers (Lumpy, Manta
and CNVkit), providing a combined output summary file and validation metrics
against NA12878 deletions. It uses the same NA12878 input as the whole genome
trio example.

To run the analysis do::

  mkdir -p NA12878-sv-eval
  cd NA12878-sv-eval
  wget https://raw.github.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-sv-getdata.sh
  bash NA12878-sv-getdata.sh
  cd work
  bcbio_nextgen.py ../config/NA12878-sv.yaml -n 16

This is large whole genome analysis and the timing and disk space requirements
for the NA12878 trio analysis above apply here as well.

RNAseq example
~~~~~~~~~~~~~~

This example aligns and creates count files for use with downstream analyses
using a subset of the SEQC data from the FDA's Sequencing Quality Control project.

Get the setup script and run it, this will download six samples from
the SEQC project, three from the HBRR panel and three from the UHRR
panel. This will require about 100GB of disk space for these input
files.  It will also set up a configuration file for the run, using
the templating system::

  wget https://raw.github.com/bcbio/bcbio-nextgen/master/config/examples/rnaseq-seqc-getdata.sh
  bash rnaseq-seqc-getdata.sh

Now go into the work directory and run the analysis::

   cd seqc/work
   bcbio_nextgen.py ../config/seqc.yaml -n 8

This will run a full scale RNAseq experiment using Tophat2 as the
aligner and will take a long time to finish on a single machine. At
the end it will output counts, Cufflinks quantitation and a set of QC
results about each lane. If you have a cluster you can `parallelize it`_
to speed it up considerably.

A nice looking standalone `report`_ of the bcbio-nextgen run can be generated using
`bcbio.rnaseq`_. Check that repository for details.

.. _templating system: https://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#automated-sample-configuration
.. _parallelize it: https://bcbio-nextgen.readthedocs.org/en/latest/contents/parallel.html
.. _bcbio.rnaseq: https://github.com/roryk/bcbio.rnaseq
.. _report: https://rawgit.com/roryk/bcbio.rnaseq/master/docs/qc-summary.html

Human genome build 38
~~~~~~~~~~~~~~~~~~~~~
Validate variant calling on human genome build 38, using two different builds
(with and without alternative alleles) and three different validation datasets
(Genome in a Bottle prepared with two methods and Illumina platinum genomes).
To run::

    mkdir -p NA12878-hg38-val
    cd NA12878-hg38-val
    wget https://raw.github.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-hg38-validate-getdata.sh
    bash NA12878-hg38-validate-getdata.sh
    mkdir -p work
    cd work
    bcbio_nextgen.py ../config/NA12878-hg38-validate.yaml -n 16

Whole genome (10x)
~~~~~~~~~~~~~~~~~~
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

- Retrieve configuration input file::

    cd config
    wget https://raw.github.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-illumina.yaml

- Run analysis on 16 core machine::

    cd work
    bcbio_nextgen.py ../config/NA12878-illumina.yaml -n 16

- Examine summary of concordance and discordance to comparison calls
  from the ``grading-summary.csv`` file in the work directory.

.. _Platinum genomes project: http://www.illumina.com/platinumgenomes/
.. _NA12878-illumina.yaml: https://raw.github.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-illumina.yaml
.. _blog post on whole genome scaling: http://bcb.io/2013/05/22/scaling-variant-detection-pipelines-for-whole-genome-sequencing-analysis/


Test suite
==========

The test suite exercises the scripts driving the analysis, so are a
good starting point to ensure correct installation. Tests use the
`pytest`_ framework. The tests are available in the bcbio source code::

     $ git clone https://github.com/bcbio/bcbio-nextgen.git

There is a small wrapper script that finds the py.test and other dependencies
pre-installed with bcbio you can use to run tests::

     $ cd tests
     $ ./run_tests.sh

You can use this to run specific test targets::

     $ ./run_tests.sh cancer
     $ ./run_tests.sh rnaseq
     $ ./run_tests.sh devel
     $ ./run_tests.sh docker

Optionally, you can run pytest directly from the bcbio install to tweak more
options. It will be in ``/path/to/bcbio/anaconda/bin/py.test``. Pass
``-s`` to ``py.test`` to see the stdout log, and ``-v`` to make py.test output more
verbose. The tests are marked with labels which you can use to run a
specific subset of the tests using the ``-m`` argument::

     $ py.test -m rnaseq

To run unit tests::

     $ py.test tests/unit

To run integration pipeline tests::

     $ py.test tests/integration

To run tests which use bcbio_vm::

     $ py.test tests/bcbio_vm

To see the test coverage, add the ``--cov=bcbio`` argument to ``py.test``.

By default the test suite will use your installed system configuration
for running tests, substituting the test genome information instead of
using full genomes. If you need a specific testing environment, copy
``tests/data/automated/post_process-sample.yaml`` to
``tests/data/automated/post_process.yaml`` to provide a test-only
configuration.

.. _pytest: http://doc.pytest.org/en/latest/
