bcbio-nextgen
-------------

A python toolkit providing best-practice pipelines for fully automated
high throughput sequencing analysis. You write a high level
configuration file specifying your inputs and analysis parameters.
This input drives a parallel pipeline that handles distributed
execution, idempotent processing restarts and safe transactional
steps. The goal is to provide a shared community resource that handles
the data processing component of sequencing analysis, providing
researchers with more time to focus on the downstream biology.

Quick start
-----------

1. Install ``bcbio-nextgen`` with all tool dependencies and data files::

         wget https://raw.github.com/chapmanb/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
         python bcbio_nextgen_install.py install_directory data_directory

producing a `system configuration file`_ referencing the installed
software and data.

2. Edit a `sample configuration file`_ to describe your samples.

3. Run analysis, distributed across 8 local cores::

         bcbio_nextgen.py bcbio_system.yaml bcbio_sample.yaml -n 8

Documentation
-------------

See the `full documentation at ReadTheDocs`_.

Pipelines
---------

Variant calling
~~~~~~~~~~~~~~~

bcbio-nextgen implements configurable best-practice pipelines for SNP
and small indel calling:

-  Sequence alignment:

   - `bowtie2`_
   - `bwa`_
   - `mosaik`_
   - `novoalign`_
   
-  Base Quality Recalibration
-  Realignment around indels
-  Variant calling:

   -  `GATK Unified Genotyper`_ (supports both GATK-lite in GATK 2.3
      and commercial restricted version in GATK 2.4+)
   -  `GATK Haplotype caller`_ (part of the commercially restricted GATK 2.4+)
   -  `FreeBayes`_
   -  `samtools mpileup`_
   -  `cortex\_var`_

-  Quality filtering, using both `GATK's Variant Quality Score
   Recalibrator`_ or hard filtering.
-  Annotation of variant effects, using `snpEff`_

It follows approaches from:

- `GATK best practice`_ guidelines for variant calling
- Marth Lab's `gkno pipelines`_

Features
--------

Distributed
~~~~~~~~~~~

The pipeline runs on single multicore machines, in compute clusters
managed by LSF or SGE using `IPython parallel`_, or on the Amazon cloud.
`This tutorial`_ describes running the pipeline on Amazon with
`CloudBioLinux`_ and `CloudMan`_.

Galaxy integration
~~~~~~~~~~~~~~~~~~

The scripts can be tightly integrated with the `Galaxy`_ web-based
analysis tool. Tracking of samples occurs via a web based LIMS system,
and processed results are uploading into Galaxy Data Libraries for
researcher access and additional analysis. See the `installation
instructions for the front end`_ and a `detailed description of the full
system`_.

.. _system configuration file: https://github.com/chapmanb/bcbio-nextgen/blob/master/config/bcbio_system.yaml
.. _sample configuration file: https://github.com/chapmanb/bcbio-nextgen/blob/master/config/bcbio_sample.yaml
.. _full documentation at ReadTheDocs: https://bcbio-nextgen.readthedocs.org
.. _GATK best practice: http://gatkforums.broadinstitute.org/discussion/1186/best-practice-variant-detection-with-the-gatk-v4-for-release-2-0
.. _GATK Unified Genotyper: http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html
.. _GATK Haplotype caller: http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html
.. _FreeBayes: https://github.com/ekg/freebayes
.. _samtools mpileup: http://samtools.sourceforge.net/mpileup.shtml
.. _cortex\_var: http://cortexassembler.sourceforge.net/index_cortex_var.html
.. _GATK's Variant Quality Score Recalibrator: http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html
.. _snpEff: http://snpeff.sourceforge.net/
.. _IPython parallel: http://ipython.org/ipython-doc/dev/index.html
.. _This tutorial: http://bcbio.wordpress.com/2011/08/19/distributed-exome-analysis-pipeline-with-cloudbiolinux-and-cloudman/
.. _CloudBioLinux: http://cloudbiolinux.org
.. _CloudMan: http://wiki.g2.bx.psu.edu/Admin/Cloud
.. _Galaxy: http://galaxy.psu.edu/
.. _installation instructions for the front end: https://bitbucket.org/galaxy/galaxy-central/wiki/LIMS/nglims
.. _detailed description of the full system: http://bcbio.wordpress.com/2011/01/11/next-generation-sequencing-information-management-and-analysis-system-for-galaxy/
.. _bwa: http://bio-bwa.sourceforge.net/
.. _bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
.. _mosaik: https://github.com/wanpinglee/MOSAIK
.. _novoalign: http://www.novocraft.com
.. _gkno pipelines: http://gkno.me/pipelines.html
