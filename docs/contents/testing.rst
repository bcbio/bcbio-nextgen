Testing
=======

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

Example pipelines
=================

We supply example input configuration files for comparison purposes
and to help in understanding the pipeline:

- An input configuration for running whole gnome variant calling with
  bwa and GATK, using Illumina's platinum genomes
  (`NA12878-illumina.yaml`_).

- Example configuration for running `ensemble variant calling`_ on
  multiple exome samples (`NA12878-ensemble.yaml`_).

We plan to utilize reference materials from NIST's `Genome in a Bottle`_
initiative to develop these into a full regression pipeline to ensure
consistency of calling between releases.

.. _NA12878-illumina.yaml: https://raw.github.com/chapmanb/bcbio-nextgen/master/config/examples/NA12878-illumina.yaml
.. _NA12878-ensemble.yaml: https://raw.github.com/chapmanb/bcbio-nextgen/master/config/examples/NA12878-ensemble.yaml
.. _ensemble variant calling: http://bcbio.wordpress.com/2013/02/06/an-automated-ensemble-method-for-combining-and-evaluating-genomic-variants-from-multiple-callers/
.. _Genome in a Bottle: http://www.genomeinabottle.org/
