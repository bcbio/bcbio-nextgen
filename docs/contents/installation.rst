Installation
------------

Automated
=========

We provide an automated script that installs 3rd party analysis tools,
required genome data, python library dependencies bundled into a
virtual environment, and produces a ready to use system configuration
file::

     wget https://raw.github.com/chapmanb/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
     python bcbio_nextgen_install.py /usr/local/share/bcbio-nextgen --tooldir=/usr/local

This installs bcbio-nextgen along with third party tools and
biological data for running human variant and RNA-seq analysis.
It places genomes, indexes and associated data files in
``/usr/local/share/bcbio-nextgen`` and tools in ``/usr/local``.
The installation is highly customizable, and you can install
additional software and data later using ``bcbio_nextgen.py upgrade``.
Run ``python bcbio_nextgen_install.py`` with no arguments to see options
for configuring the installation process. Some useful arguments are:

- ``--nosudo`` For running in environments where you lack administrator
  privileges.
- ``--isolate`` Avoid updating users ``~/.bashrc`` if installing in a
  non-standard PATH. This facilitates creation of isolated modules
  without disrupting the user's environmental setup.
- ``--nodata`` Do not install genome data.

To bootstrap installation, the machine will need to have some basic
requirements:

- Python 2.6 or 2.7, with the development libraries
  installed (the python-dev or python-devel packages).
- Compilers: gcc and g++.
- The git version control system (http://git-scm.com/).

Some steps retrieve third party tools from GitHub, which can run into
issues if you're behind a proxy or block git ports. To instruct git to
use ``https://`` globally instead of ``git://``::

    $ git config --global url.https://github.com/.insteadOf git://github.com/

The automated installer creates a fully integrated environment that
allows simultaneous updates of the framework, third party tools and
biological data. This offer the advantage over manual installation of
being able to manage and evolve a consistent analysis environment as
algorithms continue to evolve and improve. The installer is flexible
enough to handle both system integrations into standard directories
like /usr/local, as well as custom isolated installations in non-root
directories.

Upgrade
=======

We use the same automated installation process for performing upgrades
of tools, software and data in place. With a recent version of
bcbio-nextgen (0.7.0+), update with::

  bcbio_nextgen.py upgrade --tools

In addition to the installation options mentioned above, tune the
upgrade with these options:

- ``-u`` Type of upgrade to do for bcbio-nextgen code. The default is
  ``stable`` but you can also specify ``development`` to get the
  latest code from GitHub and ``skip`` to only upgrade tools and data
  without the library.

- ``--toolplus`` Specify additional categories of tools to include.
  These may require manual intervention or be data intensive. You can
  specify the argument multiple times to include multiple extra
  classes of tools. Available choices are:

  - ``protected`` Install software that requires licensing for
    commerical use. This includes the latest versions of GATK, which
    need a manual download from the GATK website. The installer
    provides full directions.
  - ``data`` Data rich supplemental tools. A good example is
    GEMINI, which provides rich annotation of variant calls
    but requires download of external data sources.

- ``--genomes`` and ``--aligners`` options add additional aligner
  indexes to download and prepare. By default we prepare a minimal
  human genome setup.

- Leave out the ``--tools`` option if you don't want to upgrade third
  party tools. If using ``--tools``, it will use the same installation
  directory as specified during installation. If you're using an older
  version that has not yet went through a successful upgrade or
  installation and saved the tool directory, you should manually
  specify ``--tooldir`` for the first upgrade. You can also pass
  ``--tooldir`` to install to a different directory.

To upgrade older bcbio-nextgen versions that don't have the ``upgrade``
command, do ``bcbio_nextgen.py -u stable`` to get the latest release
code.

On a Virtual Machine
====================
If you are looking to quickly try out bcbio-nextgen on your personal
machine before installing it on your cluster, installing bcbio-nextgen
on a virtual machine is a great way to go and is dead simple to boot,
using `Vagrant`_.

OSX
~~~
- Download and install `VirtualBox`_
- Download and install `Vagrant for OSX`_
- Get installer script::

    curl -O https://raw.github.com/chapmanb/bcbio-nextgen/master/scripts/vm/osx/vagrant_osx.sh

- Run the installer and follow the instructions::

    sh vagrant_osx.sh

.. _Vagrant for OSX: http://files.vagrantup.com/packages/7ec0ee1d00a916f80b109a298bab08e391945243/Vagrant-1.2.7.dmg
.. _VirtualBox: http://download.virtualbox.org/virtualbox/4.2.16/VirtualBox-4.2.16-86992-OSX.dmg
.. _Vagrant: http://www.vagrantup.com/

Manual process
==============

The manual process does not allow the in-place updates and management
of third party tools that the automated installer make possible. It's
a more error-prone and labor intensive process. If you find you can't
use the installer we'd love to hear why to make it more amenable to
your system.

Python code
~~~~~~~~~~~

You can install the latest release code with::

      pip install --upgrade bcbio-nextgen

Or the latest development version from GitHub::

      git clone https://github.com/chapmanb/bcbio-nextgen.git
      cd bcbio-nextgen && python setup.py build && sudo python setup.py install

This requires Python 2.7. The setup script installs
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

    fab -f cloudbiolinux/fabfile.py -H localhost install_biolinux:flavor=ngs_pipeline_minimal

You can also install them manually, adjusting locations in the
``resources`` section of your ``bcbio_system.yaml`` configuration file
as needed.

-  An aligner: we support multiple aligners, including `bwa`_,
   `novoalign`_ and `bowtie2`_
-  `Picard`_ -- BAM manipulation and processing
-  `FastQC`_ -- Generation of sequencing quality reports
-  `GATK`_ -- Variant calling and BAM preparation
-  `snpEff`_ -- Identify functional consequences of variants.

The code uses a number of Python modules, installed with the code:

-  `biopython`_
-  `pysam`_
-  `ipython`_
-  `sh`_
-  `PyYAML`_
-  `logbook`_

.. _bwa: http://bio-bwa.sourceforge.net/
.. _bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
.. _novoalign: http://www.novocraft.com
.. _Picard: http://picard.sourceforge.net/
.. _FastQC: http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/
.. _GATK: http://www.broadinstitute.org/gatk/
.. _snpEff: http://sourceforge.net/projects/snpeff/
.. _biopython: http://biopython.org
.. _pysam: http://code.google.com/p/pysam/
.. _PyYAML: http://pyyaml.org/
.. _logbook: http://packages.python.org/Logbook
.. _numpy: http://www.numpy.org/
.. _CloudBioLinux: http://cloudbiolinux.org
.. _virtualenv: http://www.virtualenv.org/en/latest/
.. _ipython: http://ipython.org/
.. _sh: http://amoffat.github.com/sh/


.. _data-requirements:

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

