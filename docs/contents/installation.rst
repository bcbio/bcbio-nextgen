Installation
------------

.. _automated-install:

Automated
=========

We provide an automated script that installs third party analysis tools,
required genome data and python library dependencies for running human variant
and RNA-seq analysis, bundled into an isolated directory or virtual environment::

     wget https://raw.github.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
     python bcbio_nextgen_install.py /usr/local/share/bcbio --tooldir=/usr/local \
       --genomes GRCh37 --aligners bwa --aligners bowtie2

bcbio should install cleanly on Linux systems. For Mac OSX, we suggest
trying `bcbio-vm <https://github.com/bcbio/bcbio-nextgen-vm>`_ which runs
bcbio on :ref:`docs-cloud` or isolates all the third party tools inside a
Docker container. bcbio-vm is still a work in progress but not all of the
dependencies bcbio uses install cleanly on OSX.

With the command line above, indexes and associated data files go in
``/usr/local/share/bcbio-nextgen`` and tools are in ``/usr/local``. If you
don't have write permissions to install into the ``/usr/local`` directories you
can install in a user directory like ``~/local`` or use ``sudo chmod`` to give
your standard user permissions. Please don't run the installer with sudo or as
the root user.  Do not use directories with ``:`` in the name, it is not
POSIX compliant and will cause installation failures.

The installation is highly customizable, and you can install
additional software and data later using ``bcbio_nextgen.py upgrade``.
Run ``python bcbio_nextgen_install.py`` with no arguments to see options
for configuring the installation process. Some useful arguments are:

- ``--isolate`` Avoid updating the user's ``~/.bashrc`` if installing in a
  non-standard PATH. This facilitates creation of isolated modules
  without disrupting the user's environmental setup. Manually edit your
  `~/.bashrc` to allow bcbio runs with::

       export PATH=/path_to_bcbio/bin:$PATH

- ``--nodata`` Do not install genome data.

The machine will need to have some basic requirements for installing and running
bcbio:

- Python 2.7, Python 3.x, or Python 2.6 plus the argparse dependency.
- Basic system setup for unpacking files: tar, gzip, unzip, bzip2, xz-utils.
- The git version control system (http://git-scm.com/)
- wget for file retrieval (https://www.gnu.org/software/wget/)

Optional tool specific requirements:

- Java 1.7, needed when running GATK < 3.6 or MuTect. This must be available in
  your path so typing ``java -version`` resolves a 1.7 version. bcbio
  distributes Java 8 as part of the anaconda installation for recent versions of
  GATK and MuTect2. You can override the Java 8 installed with bcbio by setting
  ``BCBIO_JAVA_HOME=/path/to/your/javadir`` if you have the java you want in
  ``/path/to/your/javadir/bin/java``.
- An OpenGL library, like `Mesa
  <http://mesa3d.sourceforge.net/>`_ (On Ubuntu/deb systems: ``libglu1-mesa``,
  On RedHat/rpm systems: ``mesa-libGLU-devel``). This is only required for
  cancer heterogeneity analysis with BubbleTree.
- The Pisces tumor-only variant callers requires the `Microsoft .NET runtime
  <https://www.microsoft.com/net/download/linux-package-manager/rhel/runtime-current>`_.

The `bcbio-nextgen Dockerfile
<https://github.com/bcbio/bcbio-nextgen/blob/master/Dockerfile#L5>`_ contains
the packages needed to install on bare Ubuntu systems.

The automated installer creates a fully integrated environment that allows
simultaneous updates of the framework, third party tools and biological data.
This offers the advantage over manual installation of being able to manage and
evolve a consistent analysis environment as algorithms continue to evolve and
improve. Installing this way is as isolated and self-contained as possible
without virtual machines or lightweight system containers like `Docker`_. The
:ref:`upgrade-install` section has additional documentation on including
additional genome data for supported bcbio genomes. For genome builds not
included in the defaults, see the documentation on :ref:`config-custom`.
Following installation, you should edit the pre-created system configuration
file in ``/usr/local/share/bcbio-nextgen/galaxy/bcbio_system.yaml`` to match
your local system or cluster configuration (see :ref:`tuning-cores`).

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

- ``--datatarget`` Customized installed data or download additional files not
  included by default: :ref:`datatarget-install`

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
  during installation. If you're using an older version that has not yet gone
  through a successful upgrade or installation and saved the tool directory, you
  should manually specify ``--tooldir`` for the first upgrade. You can also pass
  ``--tooldir`` to install to a different directory.

- Leave out the ``--data`` option if you don't want to get any upgrades
  of associated genome data.

- Some aligners such as STAR don't have pre-built indices due to the large file
  sizes of these. You set the number of cores to use for indexing with
  ``--cores 8``.

.. _datatarget-install:

Customizing data installation
=============================

bcbio installs associated data files for sequence processing, and you're able to
customize this to install larger files or change the defaults. Use the
``--datatarget`` flag (potentially multiple times) to customize or add new
targets.

By default, bcbio will install data files for ``variation``, ``rnaseq`` and
``smallrna`` but you can sub-select a single one of these if you don't require
other analyses. The available targets are:

- ``variation`` -- Data files required for variant calling: SNPs, indels and
  structural variants. These include files for annotation like dbSNP, associated
  files for variant filtering, coverage and annotation files.
- ``rnaseq`` -- Transcripts and indices for running RNA-seq. The transcript
  files are also used for annotating and prioritizing structural variants.
- ``smallrna`` -- Data files for doing small RNA analysis.
- ``gemini`` -- The `GEMINI <http://gemini.readthedocs.org/>`_ framework
  associates publicly available metadata with called variants, and provides
  utilities for query and analysis. This target installs the required GEMINI
  data files, including `ExAC <http://exac.broadinstitute.org/>`_.
- ``gnomad`` -- `gnomAD <http://gnomad.broadinstitute.org/>`_ is a large scale
  collection of genome variants, expanding on ExAC to include whole genome and
  more exome inputs. This is a large 25Gb download, available for human genome
  builds GRCh37, hg19 and hg38.
- ``vep`` -- Data files for the `Variant Effects Predictor (VEP)
  <http://www.ensembl.org/info/docs/tools/vep/index.html>`_. To use VEP as an
  alternative to the default installed snpEff, set ``vep`` in the
  :ref:`variant-config` configuration.
- ``dbnsfp`` Like CADD, `dbNSFP <https://sites.google.com/site/jpopgen/dbNSFP>`_
  provides integrated and generalized metrics from multiple sources to help with
  prioritizing variations for follow up. The files are large: dbNSFP is 10Gb,
  expanding to 100Gb during preparation.
- ``dbscsnv`` `dbscSNV <https://sites.google.com/site/jpopgen/dbNSFP>`_
  includes all potential human SNVs within splicing consensus regions
  (−3 to +8 at the 5’ splice site and −12 to +2 at the 3’ splice site), i.e. scSNVs,
  related functional annotations and two ensemble prediction scores for predicting their potential of altering splicing.
- ``battenberg`` Data files for `Battenberg
  <https://github.com/cancerit/cgpBattenberg>`_, which detects subclonality and
  copy number changes in whole genome cancer samples.
- ``kraken`` Database for `Kraken <https://ccb.jhu.edu/software/kraken/>`_,
  optionally used for contamination detection.
- ``ericscript`` Database for `EricScript <https://sites.google.com/site/bioericscript/>`_,
  based gene fusion detection. Supports hg38, hg19 and GRCh37.

For somatic analyses, bcbio includes `COSMIC <http://cancer.sanger.ac.uk/cosmic>`_
v68 for hg19 and GRCh37 only. Due to license restrictions, we cannot include
updated versions of this dataset and hg38 support with the installer. To prepare
these datasets yourself you can use `a utility script shipped with cloudbiolinux
<https://github.com/chapmanb/cloudbiolinux/blob/master/utils/prepare_cosmic.py>`_
that downloads, sorts and merges the VCFs, then copies into your bcbio installation::

    export COSMIC_USER="your@registered.email.edu"
    export COSMIC_PASS="cosmic_password"
    bcbio_python prepare_cosmic.py 89 /path/to/bcbio

``/path/to/bcbio/`` here is the directory one up from the ``genomes`` directory.

.. _toolplus-install:

Extra software
==============

We're not able to automatically install some useful tools due to licensing
restrictions, so we provide a mechanism to manually download and add these to
bcbio-nextgen during an upgrade with the ``--toolplus`` command line.

GATK and MuTect/MuTect2
~~~~~~~~~~~~~~~~~~~~~~~

bcbio includes an installation of GATK4, which is freely available for all uses.
This is the default runner for HaplotypeCaller or MuTect2. If you want to use an
older version of GATK, it requires manual installation. This is freely available
for academic users, but requires a `license for commercial use
<https://www.broadinstitute.org/gatk/about/#licensing>`_. It is not freely
redistributable so requires a manual download from the `GATK download`_ site.
You also need to include ``tools_off: [gatk4]`` in your configuration for runs:
see :ref:`config-changing-defaults`.

To install GATK3, register with the pre-installed gatk bioconda wrapper::

   gatk3-register /path/to/GenomeAnalysisTK.tar.bz2

If you're not using the most recent post-3.6 version of GATK, or using a nightly
build, you can add ``--noversioncheck`` to the command line to skip comparisons
to the GATK version.

`MuTect2 <https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php>`_ is distributed with GATK in versions 3.5 and later.

To install versions of GATK < 3.6, download and unzip the latest version from
the GATK distribution. Then make this jar available to bcbio-nextgen with::

    bcbio_nextgen.py upgrade --tools --toolplus gatk=/path/to/gatk/GenomeAnalysisTK.jar

This will copy the jar and update your bcbio_system.yaml and manifest files to
reflect the new version.

MuTect also has similar licensing terms and requires a license for commercial
use. After `downloading the MuTect jar
<https://www.broadinstitute.org/gatk/download/>`_, make it available to bcbio::

    bcbio_nextgen.py upgrade --tools --toolplus mutect=/path/to/mutect/mutect-1.1.7.jar

Note that muTect does not provide an easy way to query for the current version,
so your input jar needs to include the version in the name.

.. _FreeBayes and GATK comparison: http://bcb.io/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/
.. _GATK download: https://software.broadinstitute.org/gatk/download/archive


System requirements
===================

bcbio-nextgen provides a wrapper around external tools and data, so the actual
tools used drive the system requirements. For small projects, it should install
on workstations or laptops with a couple Gb of memory, and then scale as needed
on clusters or multicore machines.

Disk space requirements for the tools, including all system packages are under
4Gb. Biological data requirements will depend on the genomes and aligner indices
used, but a suggested install with GRCh37 and bowtie/bwa2 indexes uses
approximately 35Gb of storage during preparation and ~25Gb after::

    $ du -shc genomes/Hsapiens/GRCh37/*
    3.8G  bowtie2
    5.1G  bwa
    3.0G  rnaseq-2014-05-02
    3.0G  seq
    340M  snpeff
    4.2G  variation
    4.4G  vep
    23.5G total

Troubleshooting
===============

Proxy or firewall problems
~~~~~~~~~~~~~~~~~~~~~~~~~~

Some steps retrieve third party tools from GitHub, which can run into
issues if you're behind a proxy or block git ports. To instruct git to
use ``https://`` globally instead of ``git://``::

    $ git config --global url.https://github.com/.insteadOf git://github.com/

GATK or Java Errors
~~~~~~~~~~~~~~~~~~~
Most software tools used by bcbio require Java 1.8. bcbio distributes an OpenJDK
Java build and uses it so you don't need to install anything. Older versions of
GATK (< 3.6) and MuTect require a locally installed Java 1.7. If you
have version incompatibilities, you'll see errors like::

    Unsupported major.minor version 51.0

Fixing this requires either installing Java 1.7 for old GATK and MuTect or
avoiding pointing to an incorrect java (``unset JAVA_HOME``). You can also tweak
the java used by bcbio, described in the :ref:`automated-install` installation
section.

ImportErrors
~~~~~~~~~~~~
Import errors with tracebacks containing Python libraries outside of the bcbio
distribution (``/path/to/bcbio/anaconda``) are often due to other conflicting
Python installations. bcbio tries to isolate itself as much as possible but
external libraries can get included during installation due to the
PYTHONHOME or PYTHONPATH environmental variables or local site libraries.
These commands will temporary unset those to get bcbio installed, after which it
should ignore them automatically::

    $ unset PYTHONHOME
    $ unset PYTHONPATH
    $ export PYTHONNOUSERSITE=1

Finally, having a ``.pydistutils.cfg`` file in your home directory can mess with
where the libraries get installed. If you have this file in your
home directory, temporarily renaming it to something else may fix
your installation issue.

Manual process
==============

The manual process does not allow the in-place updates and management of third
party tools that the automated installer makes possible. It's a more error-prone
and labor intensive process. If you find you can't use the installer we'd love
to hear why to make it more amenable to your system. If you'd like to develop
against a bcbio installation, see the documentation on setting up a
:ref:`code-devel-infrastructure`.

Tool Requirements
~~~~~~~~~~~~~~~~~

The code drives a number of next-generation sequencing analysis tools
that you need to install on any machines involved in the processing. The
`CloudBioLinux`_ toolkit provides automated scripts to help with installation
for both software and associated data files::

    fab -f cloudbiolinux/fabfile.py -H localhost install_biolinux:flavor=ngs_pipeline_minimal

You can also install them manually, adjusting locations in the ``resources``
section of your ``bcbio_system.yaml`` configuration file as needed. The
CloudBioLinux infrastructure provides a full list of third party software
installed with bcbio-nextgen in `packages-conda.yaml`_, which lists all third
party tools installed through `Bioconda <https://bioconda.github.io/>`_

.. _CloudBioLinux: http://cloudbiolinux.org

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
