Installation
------------

You can install the latest release code with::

      pip install --upgrade bcbio-nextgen

Or the latest development version from GitHub::

      git clone https://github.com/chapmanb/bcbb.git
      cd bcbb/nextgen && python setup.py build && sudo python setup.py install

This requires either Python 2.6 or 2.7. The setup script installs
required Python library dependencies. If you'd like to install the
programs and libraries locally instead of globally, `virtualenv`_
creates an isolated, local Python installation that does not require
system install privileges.

Requirements
~~~~~~~~~~~~

The code drives a number of next-generation sequencing analysis tools
that you need to install on any machines involved in the processing. The
`CloudBioLinux`_ provides automated scripts to help with installation
for both software and associated data files. You can also install them
manually.

-  An aligner: we support multiple aligners, including `bwa`_,
   `novoalign`_ and `bowtie2`_
-  `Picard`_ -- BAM manipulation and processing
-  `FastQC`_ -- Generation of sequencing quality reports
-  `GATK`_ -- Variant calling and BAM preparation
-  `snpEff`_ -- Identify functional consequences of variants.
-  LaTeX and pdflatex for report generation

The code uses a number of Python modules, installed with the code:

-  `biopython`_
-  `numpy`_
-  `pysam`_
-  `ipython`_
-  `sh`_
-  `mako`_
-  `PyYAML`_
-  `logbook`_
-  `celery`_

.. _bwa: http://bio-bwa.sourceforge.net/
.. _bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
.. _novoalign: http://www.novocraft.com
.. _Picard: http://picard.sourceforge.net/
.. _FastQC: http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/
.. _GATK: http://www.broadinstitute.org/gatk/
.. _snpEff: http://sourceforge.net/projects/snpeff/
.. _biopython: http://biopython.org
.. _pysam: http://code.google.com/p/pysam/
.. _mako: http://www.makotemplates.org/
.. _PyYAML: http://pyyaml.org/
.. _logbook: http://packages.python.org/Logbook
.. _celery: http://celeryproject.org/
.. _numpy: http://www.numpy.org/
.. _CloudBioLinux: http://cloudbiolinux.org
.. _virtualenv: http://www.virtualenv.org/en/latest/
.. _ipython: http://ipython.org/
.. _sh: http://amoffat.github.com/sh/
