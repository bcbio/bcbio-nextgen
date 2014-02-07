Parallel execution
------------------

The pipeline runs in parallel in two different ways:

-  multiple cores -- Analyses will run in parallel using multiple cores
   on a single machine. This requires only the ``mulitprocessing``
   Python library, included by default with most Python installations.

-  parallel messaging -- This allows scaling beyond the cores
   available on a single machine, and requires multiple machines
   with a shared filesystem like standard cluster environments.
   Machine to machine communication occurs via messaging, using the
   `IPython parallel`_ framework.

Multiple cores
~~~~~~~~~~~~~~
Running using multiple cores only requires setting the ``-n``
command line flag::

    bcbio_nextgen.py bcbio_sample.yaml -t local -n 12

IPython parallel
~~~~~~~~~~~~~~~~

`IPython parallel`_ provides a distributed framework for performing
parallel computation in standard cluster environments. The
bcbio-nextgen setup script installs both IPython and `pyzmq`_, which
provides Python bindings for the `ZeroMQ`_ messaging library. The only
additional requirement is that the work directory where you run the
analysis is accessible to all processing nodes. This is typically
accomplished with a distributed file system like
`NFS`_, `Gluster`_ or `Lustre`_.

Run an analysis using ipython for parallel execution::

    bcbio_nextgen.py bcbio_sample.yaml -t ipython -n 12 -s lsf -q queue

The ``-s`` flag specifies a type of scheduler to use ``(lsf, sge, torque)``.

The ``-q`` flag specifies the queue to submit jobs to.

The ``-n`` flag defines the total number of cores to use on the
cluster during processing. The framework will select the appropriate
number of cores and type of cluster (single core versus multi-core) to
use based on the pipeline stage (see the :ref:`internals-parallel`
section in the internals documentation for more details). For
multiple core steps, the number of cores to use for programs like
``bwa``, ``novoalign`` and ``gatk`` comes from the
:ref:`config-resources` section of the configuration.
Ensure the ``cores`` specification matches the physical cores
available on machines in your cluster, and the pipeline will divide
the total cores specified by ``-n`` into the appropriate number of
multicore jobs to run.

The pipeline default parameters assume a system with minimal time to
obtain processing cores and consistent file system accessibility. These
defaults allow the system to fail fast in the case of cluster issues
which need diagnosis. For running on shared systems with high resource
usage and potential failures due to intermittent cluster issues, there
are turning parameters that increase resiliency. The ``--timeout``
flag specifies the numbers of minutes to wait for a cluster to start
up before timing out. This defaults to 15 minutes. The ``--retries``
flag specify the number of times to retry a job on failure. In systems
with transient distributed file system hiccups like lock errors or disk
availability, this will provide recoverability at the cost of
resubmitting jobs that may have failed for reproducible reasons.

Finally, the ``-r resources`` flag specifies resource options to pass along
to the underlying queue scheduler. This currently supports SGE's
``-l`` parameter, Torque's ``-l`` parameter, LSF and SLURM native flags. This allows specification
or resources to the scheduler (see the `qsub man page`_). You may specify multiple
resources, so ``-r mem=4g -r ct=01:40:00``
translates to ``-l mem=4g -l ct=01:40:00`` when passed to ``qsub`` or
``-r "account=a2010002;timelimit=04:00:00"`` when using SLURM, for
instance. SLURM and Torque support specification of an account parameter with
``-r account=your_name``, which IPython transfers into ``-A``.

Specify the `SGE parallel environment`_ to use for submitting multicore jobs
with ``-r pename=your_pe``. Since this setup
is system specific it is hard to write general code for find a
suitable environment. Specifically, when there are multiple usable
parallel environments, it will select the first one which may not be
correct. Manually specifying it with a ``pename=`` flag to resources
will ensure correct selection of the right environment.

.. _qsub man page: http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html
.. _IPython parallel: http://ipython.org/ipython-doc/dev/index.html
.. _pyzmq: https://github.com/zeromq/pyzmq
.. _ZeroMQ: http://www.zeromq.org/
.. _Gluster: http://www.gluster.org/
.. _Lustre: http://wiki.lustre.org/index.php/Main_Page
.. _NFS: https://en.wikipedia.org/wiki/Network_File_System_%28protocol%29
.. _SGE parallel environment: https://blogs.oracle.com/templedf/entry/configuring_a_new_parallel_environment

Troubleshooting
===============
Networking problems on clusters can prevent the IPython parallelization
framework from working properly. Be sure that the compute nodes on your
cluster are aware of IP addresses that they can use to communicate
with each other (usually these will be local IP addresses). Running::

    python -c 'import socket; print socket.gethostbyname(socket.gethostname())'
    
Should return such an IP address (as opposed to localhost). This can be
fixed by adding an entry to the hosts file.

The line::

    host-ip hostname
    
where ``host-ip`` is replaced by the actual IP address of the machine
and `hostname` by the machine's own hostname, should be aded to ``/etc/hosts``
on each compute node. This will probably involve contacting your local
cluster administrator.

.. _memory-management:

Memory management
~~~~~~~~~~~~~~~~~

The memory information specified in the system configuration
:ref:`config-resources` enables scheduling of memory intensive
processes. bcbio-nextgen handles memory scheduling by:

- Determining available cores and memory per machine. It uses the
  local machine for multicore runs. For parallel runs, it spawns a job
  on the schedule queue and extracts the system information from that
  machine. This expects a homogeneous set of machines within a
  cluster queue.

- Calculating the memory and core usage.
  The system configuration :ref:`config-resources` contains the
  expected core and memory usage of external programs.

- Adjusting the specified number of total cores to avoid
  over-scheduling memory. This allows running programs with more than
  the available memory per core without getting out of memory system
  errors.

- Passing total memory usage along to schedulers. The Torque, SGE and
  SLURM schedulers use this information to allocate memory to
  processes, avoiding issues with other scheduled programs using
  available memory on a shared machine.

As a result of these calculations, the cores used during processing
will not always correspond to the maximum cores provided in the input
`-n` parameter. The goal is rather to intelligently maximize cores and
memory while staying within system resources. Note that memory
specifications are for a single core, and the pipeline takes care of
adjusting this to actual cores used during processing.

Tuning systems for scale
~~~~~~~~~~~~~~~~~~~~~~~~

bcbio-nextgen scales out on clusters including hundreds of cores and is
stress tested on systems with 1000 simultaneous processes. Scaling up
often requires system specific tuning to handle simultaneous
processes. This section collects useful tips and tricks for managing
scaling issues.

Open file handles
=================

A common failure mode is having too many open file handles. This
error report can come from the IPython infrastructure logs as ZeroMQ
attempts to open sockets, or from the processing logs as third party
software gets file handles. You can check your available file handles
with ``ulimit -a | grep open``. Setting open file handle limits is
open system and cluster specific and below are tips for specific
setups.

For a Ubuntu system, edit ``/etc/security/limits.conf`` to set the
soft and hard ``nofile`` descriptors, and edit
``/etc/pam.d/common-session`` to add ``pam_limits.so``. See
`this blog post`_ for more details.

SGE needs configuration at the qmaster level. Invoke ``qconf -mconf``
from a host with admin privileges, and edit ``execd_params``::

    execd_params                 S_DESCRIPTORS=20000

.. _this blog post: https://viewsby.wordpress.com/2013/01/29/ubuntu-increase-number-of-open-files/

IO and Network File Systems
===========================

bcbio-nextgen makes use of distributed network file systems to manage
sharing large files between compute nodes. While we strive to minimize
disk-based processing by making use of pipes, the pipeline still has a
major IO component. To help manage IO and network bottlenecks, this
section contains pointers on deployments and benchmarking. Please
contribute your tips and thoughts.

- Harvard and Dell: See the 'Distributed File Systems` section of our
  `post on scaling bcbio-nextgen`_ for details about the setup within
  `Harvard FAS Research Computing`_ and thoughts on scaling and
  hardware. We also collaborate with Dell to
  test the pipeline on `Dell's Active Infrastructure for Life Sciences`_.
  We found the biggest initial factor limiting scaling was network
  bandwidth between compute and storage nodes.

.. _post on scaling bcbio-nextgen: http://bcbio.wordpress.com/2013/05/22/scaling-variant-detection-pipelines-for-whole-genome-sequencing-analysis/
.. _Harvard FAS Research Computing: http://rc.fas.harvard.edu/
.. _Dell's Active Infrastructure for Life Sciences: http://dell.com/ai-hpc-lifesciences

Cloud support
~~~~~~~~~~~~~

`Amazon Web Services`_ provide a flexible cloud based environment for
running analyses. Cloud approaches offer the ability to perform
analyses at scale with no investment in local hardware. In addition to
the potential advantages for traditional cluster users, shared images
on top of this infrastructure can make these analysis pipelines
available to anyone. `This tutorial`_ describes running the pipeline
on Amazon with `CloudBioLinux`_ and `CloudMan`_.

The scripts can also be tightly integrated with the `Galaxy`_ web-based
analysis tool. Tracking of samples occurs via a web based LIMS system,
and processed results are uploading into Galaxy Data Libraries for
researcher access and additional analysis. See the `installation
instructions for the front end`_ and a `detailed description of the full
system`_.

.. _Amazon Web Services: https://aws.amazon.com/
.. _This tutorial: http://bcbio.wordpress.com/2011/08/19/distributed-exome-analysis-pipeline-with-cloudbiolinux-and-cloudman/
.. _CloudBioLinux: http://cloudbiolinux.org
.. _CloudMan: http://wiki.g2.bx.psu.edu/Admin/Cloud

.. _Galaxy: http://galaxy.psu.edu/
.. _installation instructions for the front end: https://bitbucket.org/galaxy/galaxy-central/wiki/LIMS/nglims
.. _detailed description of the full system: http://bcbio.wordpress.com/2011/01/11/next-generation-sequencing-information-management-and-analysis-system-for-galaxy/
