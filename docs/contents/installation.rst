Requirements
------------

The code drives a number of next-generation sequencing analysis tools
that you need to install on any machines involved in the processing. The
`CloudBioLinux`_ provides automated scripts to help with installation
for both software and associated data files. You can also install them
manually.

Sequencing analysis
~~~~~~~~~~~~~~~~~~~

-  An aligner: we support multiple aligners, including `bwa`_,
   `novoalign`_ and `bowtie2`_
-  `Picard`_
-  `FastQC`_
-  `GATK`_
-  `snpEff`_

Processing infrastructure
~~~~~~~~~~~~~~~~~~~~~~~~~

-  RabbitMQ for communication between machines
-  LaTeX and pdflatex for report generation

Python modules installed with the package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `biopython`_
-  `pysam`_
-  `mako`_
-  `PyYAML`_
-  `logbook`_
-  `celery`_
-  `numpy`_

Installation
------------

Scripts and configuration
~~~~~~~~~~~~~~~~~~~~~~~~~

Clone a copy of the code from GitHub:

::

      git clone https://github.com/chapmanb/bcbb.git

Use a recent version of Python 2 (2.6 or 2.7), and install with:

::

      cd bcbb/nextgen && python setup.py install

The setup script installs required python library dependencies.

Copy the YAML & ini files in config and adjust them to match your
environment. It is also a good idea to set your $PATH pointing to any
third-party binaries you are using.

Parallel execution
~~~~~~~~~~~~~~~~~~

The pipeline runs in parallel in two different ways:

-  multiple cores -- Analyses will run in parallel using multiple cores
   on a single machine. This requires only the ``mulitprocessing``
   Python library, included by default with most Python installations.
   Change ``num_cores`` in your ``post_process.yaml`` file to specify
   the number of parallel cores to use.

-  parallel messaging -- This allows scaling beyond the cores on a
   single machine. It requires multiple machines with a shared
   filesystem, with communication handled using RabbitMQ messaging. This
   is ideally suited for clusters.

To enable parallel messaging:

1. Configure RabbitMQ as described below. Ensure all processing machines
   can talk to the RabbitMQ server on port 5672. Update
   ``universe_wsgi.ini`` to contain the server details.

2. Edit your ``post_process.yaml`` file to set parameters in the
   ``distributed`` section corresponding to your environment: this
   includes the type of cluster management and arguments to start jobs.

3. Run ``bcbio_nextgen.py`` with parameters for a distributed cluster
   environment. It takes care of starting worker nodes, running the
   processing, and then cleaning up after jobs:

   bcbio\_nextgen.py post\_process.yaml flowcell\_dir run\_info.yaml -t
   messaging -n 20

Virtual development environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The installation instructions assume that you have full root access to
install python modules and packages (production environment). If this is
not the case, you may want to install a python VirtualEnv and other
tools automatically on your $HOME to ease your development needs using
the following script:

::

        wget https://bitbucket.org/brainstorm/custom_env/raw/1cd4f4ae27d5/pyHost.sh && ./pyHost.sh

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
