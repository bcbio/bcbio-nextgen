Installation
------------

We provide an automated script that installs 3rd party analysis tools,
required genome data, python library dependencies bundled into a
virtual environment, and produces a ready to use system configuration
file::

     wget https://raw.github.com/chapmanb/bcbb/master/nextgen/scripts/bcbio_nextgen_install.py
     python bcbio_nextgen_install.py install_directory data_directory

By default the script downloads genomes, indexes and associated data
files for human variant and RNA-seq analysis. Run
``python bcbio_nextgen_install.py`` with no arguments to see options
for configuring the installation process. There is a ``--nosudo``
argument for running in environments where you lack administrator
privileges. 

This script requires that you can do a ``ssh localhost`` to your
installation machine. If you'd like to do this without any passwords
you can setup your ssh keys with::

    $ cat ~/.ssh/id_rsa.pub >> ~/.ssh/authorized_keys

If you'd prefer more control over installation, follow the manual
steps for installing each component detailed below.

Python code
~~~~~~~~~~~

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

Tool Requirements
~~~~~~~~~~~~~~~~~

The code drives a number of next-generation sequencing analysis tools
that you need to install on any machines involved in the processing. The
`CloudBioLinux`_ toolkit provides automated scripts to help with installation
for both software and associated data files::

    fab -f cloudbiolinux/fabfile.py -H localhost install_biolinux:flavor=ngs_pipeline

You can also install them manually, adjusting locations in the
``resources`` section of your ``bcbio_system.yaml`` configuration file
as needed.

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

Data requirements
~~~~~~~~~~~~~~~~~

In addition to existing bioinformatics software the pipeline requires
associated data files for reference genomes, including pre-built indexes
for aligners. The `CloudBioLinux`_ toolkit again provides an automated
way to download and prepare these reference genomes::

    fab -f data_fabfile.py -H localhost -c your_fabricrc.txt install_data_s3:your_biodata.yaml

The `biodata.yaml`_ file contains information about what genomes to
download. The `fabricrc.txt`_ describes where to install the genomes
by adjusting the ``data_files`` variable. This creates a tree
structure that includes a set of Galaxy-style location files to
describe locations of indexes::

    ├── galaxy
    │   ├── tool-data
    │   │   ├── alignseq.loc
    │   │   ├── bowtie_indices.loc
    │   │   ├── bwa_index.loc
    │   │   ├── sam_fa_indices.loc
    │   │   └── twobit.loc
    │   └── tool_data_table_conf.xml
    ├── genomes
    │   ├── Hsapiens
    │   │   ├── GRCh37
    │   │   └── hg19
    │   └── phiX174
    │       └── phix
    └── liftOver

Individual genome directories contain indexes for aligners in
individual sub-directories prefixed by the aligner name. This
structured scheme helps manage aligners that don't have native Galaxy
`.loc` files. The automated installer will download and set this up
automatically::

    `-- phix
        |-- bowtie
        |   |-- phix.1.ebwt
        |   |-- phix.2.ebwt
        |   |-- phix.3.ebwt
        |   |-- phix.4.ebwt
        |   |-- phix.rev.1.ebwt
        |   `-- phix.rev.2.ebwt
        |-- bowtie2
        |   |-- phix.1.bt2
        |   |-- phix.2.bt2
        |   |-- phix.3.bt2
        |   |-- phix.4.bt2
        |   |-- phix.rev.1.bt2
        |   `-- phix.rev.2.bt2
        |-- bwa
        |   |-- phix.fa.amb
        |   |-- phix.fa.ann
        |   |-- phix.fa.bwt
        |   |-- phix.fa.pac
        |   |-- phix.fa.rbwt
        |   |-- phix.fa.rpac
        |   |-- phix.fa.rsa
        |   `-- phix.fa.sa
        |-- novoalign
        |   `-- phix
        |-- seq
        |   |-- phix.dict
        |   |-- phix.fa
        |   `-- phix.fa.fai
        `-- ucsc
            `-- phix.2bit
    
.. _fabricrc.txt: https://github.com/chapmanb/cloudbiolinux/blob/master/config/fabricrc.txt
.. _biodata.yaml: https://github.com/chapmanb/cloudbiolinux/blob/master/config/biodata.yaml
    
