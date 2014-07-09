Presentations
=============

Bioinformatics Open Source Conference 2013
------------------------------------------

-  `BOSC 2013 Slides`_
-  `BOSC 2013 Video`_
-  `BOSC 2013 Conference website`_

Arvados Summit 2013
-------------------

-  `Arvados Summit Slides`_
-  `Arvados Summit website`_

Scientific Python 2013
----------------------

-  `SciPy 2013 Video`_
-  `SciPy 2013 Conference website`_

Genome Informatics 2013
-----------------------

-  `GI 2013 Presentation slides`_

Mt Sinai: Strategies for accelerating the genomic sequencing pipeline
---------------------------------------------------------------------

- `Mt Sinai workshop slides`_
- `Mt Sinai workshop website`_

Non-conference talks
--------------------

- `Novartis slides`_ (21 January 2014)

Abstract
~~~~~~~~

**Community Development of Validated Variant Calling Pipelines**

*Brad Chapman, Rory Kirchner, Oliver Hofmann and Winston Hide Harvard
School of Public Health, Bioinformatics Core, Boston, MA, 02115*

Translational research relies on accurate identification of genomic
variants. However, rapidly changing best practice approaches in
alignment and variant calling, coupled with large data sizes, make it a
challenge to create reliable and reproducible variant calls. Coordinated
community development can help overcome these challenges by sharing
testing and updates across multiple groups. We describe bcbio-nextgen, a
distributed multi-architecture pipeline that automates variant calling,
validation and organization of results for query and visualization. It
creates an easily installable, reliable infrastructure from
best-practice open source tools with the following goals:

-  **Quantifiable:** Validates variant calls against known reference
   materials developed by the `Genome in a Bottle`_ consortium. The
   `bcbio.variation`_ toolkit automates scoring and assessment of calls
   to identify regressions in variant identification as calling
   pipelines evolve. Incorporation of multiple variant calling
   approaches from `Broad's GATK best practices`_ and the `Marth lab's
   gkno software`_ enables informed comparisons between current and
   future algorithms.

-  **Scalable:** bcbio-nextgen handles large population studies with
   hundreds of whole genome samples by parallelizing on a wide variety
   of schedulers and multicore machines, setting up different ad hoc
   cluster configurations for each workflow step. Work in progress
   includes integration with virtual environments, including `Amazon Web
   Services`_ and `OpenStack`_.

-  **Accessible:** Results automatically feed into tools for query and
   investigation of variants. The `GEMINI framework`_ provides a
   queryable database associating variants with a wide variety of genome
   annotations. The `o8`_ web-based tool visualizes the work of variant
   prioritization and assessment.

-  **Community developed:** bcbio-nextgen is widely used in multiple
   sequencing centers and research laboratories. We actively encourage
   contributors to the code base and make it easy to get started with a
   fully automated installer and updater that prepares all third party
   software and reference genomes.

Links from the presentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `HugeSeq`_
-  `Genome Comparison & Analytic Testing`_ at Bioplanet
-  `Peter Block’s “Community” book`_
-  `CloudBioLinux`_ and `Homebrew Science`_ as installation frameworks;
   `Conda`_ as Python environment
-  bcbio `documentation`_ at ReadTheDocs
-  `Arvados framework`_ for meta data tracking, NGS processing and data
   provenance
-  Notes on `improved scaling for NGS workflows`_
-  Genomic Reference Materials from `Genome in a Bottle`_
-  Comparison of `aligners and callers`_ using NIST reference materials
-  Callers and `minimal BAM preparation workflows`_
-  `Coverage assessment`_

.. _BOSC 2013 Slides: http://chapmanb.github.io/bcbb/talks/bosc2013_bcbio_nextgen/chapmanb_bosc2013_bcbio.html#/
.. _BOSC 2013 Video: http://www.youtube.com/watch?v=dT5UEU0xF1Q
.. _BOSC 2013 Conference website: http://www.open-bio.org/wiki/BOSC_2013
.. _Arvados Summit Slides: https://github.com/chapmanb/bcbb/raw/master/talks/arvados2013_bcbio_nextgen/chapman_arvadossum_bcbio.pdf
.. _Arvados Summit website: https://arvados.org/projects/arvados/wiki/Arvados_Summit_-_Fall_2013
.. _SciPy 2013 Video: https://www.youtube.com/watch?v=qNMPh0pIpBE
.. _SciPy 2013 Conference website: https://conference.scipy.org/scipy2013/
.. _GI 2013 Presentation slides: https://dl.dropboxusercontent.com/u/407047/Work/Presentations/20131102%20CSHL%20Genome%20Informatics/20131101%20CSHL%20GI2013%20bcbio.pdf
.. _Genome in a Bottle: http://www.genomeinabottle.org/
.. _bcbio.variation: https://github.com/chapmanb/bcbio.variation
.. _Broad's GATK best practices: http://gatkforums.broadinstitute.org/discussion/1186/best-practice-variant-detection-with-the-gatk-v4-for-release-2-0
.. _Marth lab's gkno software: http://gkno.me/
.. _Amazon Web Services: https://aws.amazon.com/
.. _OpenStack: http://www.openstack.org/
.. _GEMINI framework: https://github.com/arq5x/gemini#readme
.. _o8: https://github.com/chapmanb/o8#readme
.. _HugeSeq: http://github.com/StanfordBioinformatics/HugeSeq
.. _Genome Comparison & Analytic Testing: http://www.bioplanet.com/gcat
.. _Peter Block’s “Community” book: http://www.amazon.com/Community-Structure-Belonging-Peter-Block/dp/1605092770
.. _CloudBioLinux: http://cloudbiolinux.org/
.. _Homebrew Science: https://github.com/Homebrew/homebrew-science
.. _Conda: http://www.continuum.io/blog/conda
.. _documentation: bcbio-nextgen.readthedocs.org
.. _Arvados framework: https://arvados.org/
.. _improved scaling for NGS workflows: http://bcbio.wordpress.com/2013/05/22/scaling-variant-detection-pipelines-for-whole-genome-sequencing-analysis/
.. _aligners and callers: http://bcbio.wordpress.com/2013/05/06/framework-for-evaluating-variant-detection-methods-comparison-of-aligners-and-callers/
.. _minimal BAM preparation workflows: http://bcbio.wordpress.com/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/
.. _Coverage assessment: https://github.com/chapmanb/bcbio.coverage
.. _Mt Sinai workshop website: http://www.hpcwire.com/event/strategies-accelerating-genomic-sequencing-pipeline/
.. _Mt Sinai workshop slides: https://github.com/chapmanb/bcbb/raw/master/talks/mtsinai2013_bcbio_nextgen/chapman_mtsinai_bcbio.pdf
.. _Novartis slides: https://github.com/chapmanb/bcbb/raw/master/talks/novartis2014_bcbio_nextgen/chapman_bcbio.pdf
