Installation
------------

Automated
=========

We provide an automated script that installs third party analysis tools,
required genome data and python library dependencies for running human variant
and RNA-seq analysis, bundled into an isolated directory or virtual environment::

     wget https://raw.github.com/chapmanb/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
     python bcbio_nextgen_install.py /usr/local/share/bcbio --tooldir=/usr/local \
       --genomes GRCh37 --aligners bwa --aligners bowtie2

It places genomes, indexes and associated data files in
``/usr/local/share/bcbio-nextgen`` and tools in ``/usr/local``. You should edit
the pre-created system configuration file in
``/usr/local/share/bcbio-nextgen/galaxy/bcbio_system.yaml``
to match your local system or cluster configuration.

The installation is highly customizable, and you can install
additional software and data later using ``bcbio_nextgen.py upgrade``.
Run ``python bcbio_nextgen_install.py`` with no arguments to see options
for configuring the installation process. Some useful arguments are:

- ``--sudo`` Enable installation in privileged directories and allow the
  installer to update system packages.
- ``--isolate`` Avoid updating the user's ``~/.bashrc`` if installing in a
  non-standard PATH. This facilitates creation of isolated modules
  without disrupting the user's environmental setup.
- ``--nodata`` Do not install genome data.

To bootstrap installation, the machine will need to have some basic
requirements:

- Python 2.6 or 2.7, with the development libraries
  installed (the python-dev or python-devel packages).
- Compilers: Recent versions of gcc and g++. gcc 4.8.x is well tested,
  although other versions should work fine.
- The git version control system (http://git-scm.com/).
- wget for file retrieval (https://www.gnu.org/software/wget/)
- unzip
- zlib (with development libraries)

If you're not using the ``--sudo`` option, please see :ref:`isolated-install`
for additional system requirements needed to bootstrap the full system on
minimal machines. The
`bcbio-nextgen Dockerfile <https://github.com/chapmanb/bcbio-nextgen/blob/master/Dockerfile#L5>`_
contains bootstrap package information to install on bare Ubuntu systems.

The automated installer creates a fully integrated environment that
allows simultaneous updates of the framework, third party tools and
biological data. This offers the advantage over manual installation of
being able to manage and evolve a consistent analysis environment as
algorithms continue to evolve and improve. The installer is flexible
enough to handle both system integrations into standard directories
like /usr/local, as well as custom isolated installations in non-root
directories. The :ref:`upgrade-install` section has additional
documentation on including additional genome data, and the section on
:ref:`toolplus-install` describes how to add commercially restricted software
like GATK.

.. _isolated-install:

Isolated installations
======================

To install bcbio-nextgen in an isolated non-root environment::

    python bcbio_nextgen_install.py /path_to_bcbio --tooldir=/path_to_bcbio --isolate \
       --genomes GRCh37 --aligners bwa --aligners bowtie2

This requires the following additional system requirements to be in place:

- Java 1.7
- Ruby
- R with Rscript (currently optional, but increasingly used in the pipeline)
- bzip2 (with development libraries)
- curl (with development libraries)
- curses (with development libraries)

Installing this way is as isolated and self-contained as possible
without virtual machines or lightweight system containers. To ensure
access to the executables, system libraries and Perl libraries update
your `~/.bashrc` with::

    export PATH=/path_to_bcbio/bin:$PATH
    export LD_LIBRARY_PATH=/path_to_bcbio/lib:$LD_LIBRARY_PATH
    export PERL5LIB=/path_to_bcbio/lib/perl5:/path_to_bcbio/perl5/site_perl:${PERL5LIB}

This installation process is not easily re-locatable due to absolute
filesystem pointers within the installation directory. We plan to move
towards utilizing `Docker`_ containers to provide a fully isolated software
installation.

.. _Docker: http://www.docker.io/

.. _upgrade-install:

Upgrade
=======

We use the same automated installation process for performing upgrades
of tools, software and data in place. Since there are multiple targets
and we want to avoid upgrading anything unexpectedly, we have specific
arguments for each. Generally, you'd want to upgrade the code, tools
and data together with::

  bcbio_nextgen.py upgrade -u stable --tools --data

Tune the upgrade with these options:

- ``-u`` Type of upgrade to do for bcbio-nextgen code. ``stable``
  gets the most recent released version and ``development``
  retrieves the latest code from GitHub.

- ``--toolplus`` Specify additional tools to include. See the section on
  :ref:`toolplus-install` for more details.

- ``--genomes`` and ``--aligners`` options add additional aligner
  indexes to download and prepare. ``bcbio_nextgen.py upgrade -h`` lists
  available genomes and aligners. If you want to install multiple genomes or
  aligners at once, specify ``--genomes`` or ``--aligners``
  multiple times, like this:
  ``--genomes GRCh37 --genomes mm10 --aligners bwa --aligners bowtie2``

- Leave out the ``--tools`` option if you don't want to upgrade third party
  tools. If using ``--tools``, it will use the same directory as specified
  during installation. If you're using an older version that has not yet went
  through a successful upgrade or installation and saved the tool directory, you
  should manually specify ``--tooldir`` for the first upgrade. You can also pass
  ``--tooldir`` to install to a different directory.

- Leave out the ``--data`` option if you don't want to get any upgrades
  of associated genome data.

.. _toolplus-install:

System requirements
===================

bcbio-nextgen provides a wrapper around external tools and data, so the actual
tools used drive the system requirements. For small projects, it should install
on workstations or laptops with a couple Gb of memory, and then scale as needed
on clusters or multicore machines.

Disk space requirements for the tools, including all system packages are under
4Gb. Biological data requirements will depend on the genomes and aligner indices
used, but a suggested install with GRCh37 and bowtie/bwa2 indexes uses
under 25Gb of storage::

    $ du -shc genomes/Hsapiens/GRCh37/*
    3.8G  bowtie2
    5.1G  bwa
    3.0G  rnaseq-2014-05-02
    3.0G  seq
    340M  snpeff
    4.2G  variation
    4.4G  vep
    23.5G total

.. _extra-install:

Extra software and data
=======================

We're not able to automatically install some useful tools due to licensing
restrictions, so provide a mechanism to manually download and add these to
bcbio-nextgen during an upgrade with the ``--toolplus`` command line. This also
includes mechanisms to add in large annotation files not included by default.

GATK and muTect
~~~~~~~~~~~~~~~

Calling variants with GATK's HaplotypeCaller or UnifiedGenotyper requires manual
installation of the latest GATK release. This is freely available for academic
users, but requires a manual download from the `GATK download`_ site.  Appistry
provides `a distribution of GATK for commercial users`_.  We distribute the last
freely available GATK-lite release (2.3.9) as part of the automated install but
don't recommend using this for calling. If you don't want to use the restricted
GATK version, freely available callers like FreeBayes provide a better
alternative than older GATK versions. See the `FreeBayes and GATK comparison`_
for a full evaluation.

To install GATK, download and unzip the latest version from the GATK or Appistry
distributions. Then make this jar available to bcbio-nextgen with::

    bcbio_nextgen.py upgrade --tools --toolplus gatk=/path/to/gatk/GenomeAnalysisTK.jar

This will copy the jar and update your bcbio_system.yaml and manifest files to
reflect the new version.

For muTect, we provide the latest 1.1.5 jar, but commercial users need to obtain
the Appistry muTect distribution. To make this jar available to bcbio-nextgen::

    bcbio_nextgen.py upgrade --tools --toolplus mutect=/path/to/appistry/muTect-1.1.5.jar

Note that muTect does not provide an easy way to query for the current version,
so your input jar needs to include the version in the name.

GEMINI
~~~~~~

``-- toolplus`` is also used to install data rich supplemental software which is
not installed by default such as GEMINI. We're making changes to automatically
include these tools in the default install, but for now include  GEMINI with::

    bcbio_nextgen.py upgrade --tools --toolplus data

dbNSFP and CADD
~~~~~~~~~~~~~~~

Two useful databases for evaluating the potential impact of variations are
`CADD`_ and `dbNSFP`_. They provide integrated and generalized metrics from
multiple sources to help with prioritizing variations for follow up. The files
are large: dbNSFP is 10Gb, expanding to 100Gb during preparation; and CADD is
30Gb. As a result they are not included in an install by default. You can add them,
either together or individually, using ``--toolplus``::

    bcbio_nextgen.py upgrade --tools --toolplus cadd --toolplus dbnsfp --data

When installed, GEMINI will automatically include `CADD`_ annotations as part of
the created SQLite database. Setting `VEP`_ in the :ref:`variant-config`
configuration will include annotation of VCFs with `dbNSFP`_.

Both tools are freely available for non-commercial research, but require licensing
for commercial usage.

.. _CADD: http://cadd.gs.washington.edu/home
.. _dbNSFP: https://sites.google.com/site/jpopgen/dbNSFP
.. _VEP: http://www.ensembl.org/info/docs/tools/vep/index.html
.. _GATK download: http://www.broadinstitute.org/gatk/download
.. _a distribution of GATK for commercial users: http://www.appistry.com/gatk
.. _FreeBayes and GATK comparison: http://bcbio.wordpress.com/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/

Troubleshooting
===============

Old bcbio version support
~~~~~~~~~~~~~~~~~~~~~~~~~

The upgrade approach changed slightly as of 0.7.5 to be more
consistent.  In earlier versions, to get a full upgrade leave out the
``--data`` argument since that was the default. The best approach if
you find the arguments are out of date is to do a ``bcbio_nextgen.py
upgrade -u stable`` to get the latest version, then proceed
again. Pre 0.7.0 versions won't have the ``upgrade`` command and need
``bcbio_nextgen.py -u stable`` to get up to date.

Proxy or firewall problems
~~~~~~~~~~~~~~~~~~~~~~~~~~

Some steps retrieve third party tools from GitHub, which can run into
issues if you're behind a proxy or block git ports. To instruct git to
use ``https://`` globally instead of ``git://``::

    $ git config --global url.https://github.com/.insteadOf git://github.com/


ImportError: No module named conda.cli
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Having a PYTHONHOME or PYTHONPATH set can cause installation troubles,
if you are seeing an error like the above, unsetting these two environment
variables will help. Fix that with::

    $ unset PYTHONHOME
    $ unset PYTHONPATH

Other import errors
~~~~~~~~~~~~~~~~~~~
Having a .pydistutils.cfg file in your home directory can mess with
where the libraries get installed. If you have this file in your
home directory, temporarily renaming it to something else may fix
your installation issue.

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
- Ensure your system has wget installed.
- Get and run the installer script::

    mkdir bcbio && cd bcbio
    wget https://raw.github.com/chapmanb/bcbio-nextgen/master/scripts/vm/osx/vagrant_osx.sh
    sh vagrant_osx.sh

.. _Vagrant for OSX: http://www.vagrantup.com/downloads.html
.. _VirtualBox: https://www.virtualbox.org/wiki/Downloads
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

.. _virtualenv: http://www.virtualenv.org/en/latest/

Tool Requirements
~~~~~~~~~~~~~~~~~

The code drives a number of next-generation sequencing analysis tools
that you need to install on any machines involved in the processing. The
`CloudBioLinux`_ toolkit provides automated scripts to help with installation
for both software and associated data files::

    fab -f cloudbiolinux/fabfile.py -H localhost install_biolinux:flavor=ngs_pipeline_minimal

You can also install them manually, adjusting locations in the
``resources`` section of your ``bcbio_system.yaml`` configuration file
as needed.  The CloudBioLinux infrastructure provides a full list of third party
software installed with bcbio-nextgen:

- `packages-homebrew.yaml`_ -- All third party tools installed through the
  Homebrew/Linuxbrew package manager.
- `custom.yaml`_ -- All third party tools installed via CloudBioLinux's custom
  installation procedure.

.. _CloudBioLinux: http://cloudbiolinux.org
.. _packages-homebrew.yaml: https://github.com/chapmanb/cloudbiolinux/blob/master/contrib/flavor/ngs_pipeline_minimal/packages-homebrew.yaml
.. _custom.yaml : https://github.com/chapmanb/cloudbiolinux/blob/master/contrib/flavor/ngs_pipeline_minimal/custom.yaml

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
