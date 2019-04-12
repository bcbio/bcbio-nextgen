Parallel execution
------------------

The pipeline runs in parallel in two different ways:

-  multiple cores -- Analyses will run in parallel using multiple cores
   on a single machine. This requires only the ``multiprocessing``
   Python library, included by default with most Python installations.

-  parallel messaging -- This allows scaling beyond the cores
   available on a single machine, and requires multiple machines
   with a shared filesystem like standard cluster environments.
   Machine to machine communication occurs via messaging, using the
   `IPython parallel`_ framework.

.. _tuning-cores:

Tuning core and memory usage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bcbio has two ways to specify core usage, helping provide options for
parallelizing different types of processes:

- Total available cores: specified with ``-n`` on the commandline, this tells
  bcbio how many total cores to use. This applies either to a local multicore
  run or a distributed job.

- Maximum cores to use for multicore processing of individual jobs. You specify
  this in the ``resource`` section of either a sample YAML file
  (:ref:`sample-resources`) or ``bcbio_system.yaml``. Ideally you specify this
  in the ``default`` section (along with memory usage). For example, this would
  specify that processes using multiple cores can get up to 16 cores with 2G of
  memory per core::

      resources:
        default:
          memory: 2G
          cores: 16
          jvm_opts: ["-Xms750m", "-Xmx2000m"]

bcbio uses these settings, along with memory requests, to determine how to
partition jobs. For example, if you had ``-n 32`` and ``cores: 16`` for a run on
a single 32 core machine, this would run two simultaneous bwa mapping jobs using
16 cores each.

Memory specifications (both in ``memory`` and ``jvm_opts``) are per-core. bcbio
takes care of adjusting this memory to match the cores used. In the example
above, if bcbio was running a 16 core java process, it would use 32Gb of memory
for the JVM, adjusting ``Xmx`` and ``Xms`` to match cores used. Internally bcbio
looks at the memory and CPU usage on a machine and matches your configuration
options to the available system resources. It will scale down core requests if
memory is limiting, avoiding over-scheduling resources during the run. You
ideally want to set both ``memory`` and ``jvm_opts`` to match the average memory
per core on the run machine and adjust upwards if this does not provide enough
memory for some processes during the run.

For single machine runs with a small number of samples, you generally want to
set ``cores`` close to or equal the number of total cores you're allocating to
the job with ``-n``. This will allow individual samples to process as fast as
possible and take advantage of multicore software.

For distributed jobs, you want to set ``cores`` to match the available cores on
a single node in your cluster, then use ``-n`` as a multiple of this to
determine how many nodes to spin up. For example, ``cores: 16`` and ``-n 64``
would try to make four 16 core machines available for analysis.

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

The ``-s`` flag specifies a type of scheduler to use ``(lsf, sge, torque, slurm, pbspro)``.

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
``-l`` parameter, Torque's ``-l`` parameter and LSF and SLURM native flags. This allows specification
or resources to the scheduler (see the `qsub man page`_). You may specify multiple
resources, so ``-r mem=4g -r ct=01:40:00``
translates to ``-l mem=4g -l ct=01:40:00`` when passed to ``qsub`` or
``-r "account=a2010002" -r "timelimit=04:00:00"`` when using SLURM, for
instance. SLURM and Torque support specification of an account parameter with
``-r account=your_name``, which IPython transfers into ``-A``.

SGE supports special parameters passed using resources to help handle the
heterogeneity of possible setups.

Specify an `SGE parallel environment
<https://docs.oracle.com/cd/E19957-01/820-0698/6ncdvjcmd/index.html>`_
that supports using multiple cores on a
single node with ``-r pename=your_pe``. Since this setup is system specific it
is hard to write general code for find a suitable environment. Specifically,
when there are multiple usable parallel environments, it will select the first
one which may not be correct. Manually specifying it with a ``pename=`` flag to
resources will ensure correct selection of the right environment. If you're
administering a grid engine cluster and not sure how to set this up you'd
typically want a ``smp`` queue using ``allocation_rule: $pe_slots`` like in this
`example pename configuration
<https://github.com/WGLab/biocluster/blob/431a05f6dfd532205aacfc7477ac740b0e7b2a0a/03%20System%20customization.md#setting-up-parallel-environment>`_
or `smp template <https://gist.github.com/dan-blanchard/6586533#file-smp_template>`_.

SGE has other specific flags you may want to tune, depending on your setup. To
specify an advanced reservation with the ``-ar`` flag, use ``-r ar=ar_id``. To
specify an alternative memory management model instead of ``mem_free`` use ``-r
memtype=approach``. It is further recommended to configure ``mem_free`` (or any
other chosen memory management model) as a consumable, requestable resource in
SGE to prevent overfilling hosts that do not have sufficient memory per slot.
This can be done in two steps. First, launch ``qmon`` as an admin, select
``Complex Configuration`` in qmon, click on ``mem_free`, under the
``Consumable`` dialog select ``JOB`` (instead of ``YES`` or ``NO``) and finally
click ``Modify`` for the changes to take effect. Secondly, for each host in the
queue, configure ``mem_free`` as a complex value. If a host called ``myngshost``
has 128GB of RAM, the corresponding command would be ``qconf -mattr exechost
complex_values mem_free=128G myngshost``

There are also special ``-r`` resources parameters to support pipeline configuration:

- ``-r conmem=4`` -- Specify the memory for the controller process, in Gb. This
  currently applies to SLURM processing and defaults to 4Gb.

- ``-r minconcores=2`` -- The minimum number of cores to use for the controller
  process. The controller one works on a single core but this can help in
  queues where you can only specify multicore jobs.

- ``-r mincores=16`` -- Specify the minimum number of cores to batch together
  for parallel single core processes like variant calling. This will run
  multiple processes together under a single submission to allow sharing of
  resources like memory, which is helpful when a small percentage of the time a
  process like variant calling will use a lot of memory. By default, bcbio will
  calculate ``mincores`` based on specifications for multicore calling so this
  doesn't normally require a user to set.

.. _qsub man page: http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html
.. _IPython parallel: http://ipython.org/ipython-doc/dev/index.html
.. _pyzmq: https://github.com/zeromq/pyzmq
.. _ZeroMQ: http://www.zeromq.org/
.. _Gluster: http://www.gluster.org/
.. _Lustre: http://wiki.lustre.org/index.php/Main_Page
.. _NFS: https://en.wikipedia.org/wiki/Network_File_System_%28protocol%29

Troubleshooting
~~~~~~~~~~~~~~~
Diagnosing job failures
=======================

Parallel jobs can often terminate with rather generic failures like any of the
following:

- ``joblib/parallel.py", ... TypeError: init() takes at least 3 arguments (2 given)``
- ``Multiprocessing exception:``
- ``CalledProcessError: Command '<command line that failed>``

These errors unfortunately don't help diagnose the problem, and you'll likely
see the actual error triggering this generic exception earlier in the run. This
error can often be hard to find due to parallelization.

If you run into a confusing failure like this, the best approach is to re-run
with a single core::

    bcbio_nextgen.py your_input.yaml -n 1

which should produce a more helpful debug message right above the failure.

It's also worth re-trying the failed command line outside of bcbio to look for
errors. You can find the failing command by cross-referencing the error message
with command lines in ``log/bcbio-nextgen-commands.log``. You may have to change
temporary directories (``tx/tmp**``) in some of the job outputs. Reproducing the
error outside of bcbio is a good first step to diagnosing and fixing the
underlying issue.

No parallelization where expected
=================================

This may occur if the current execution is a re-run of a previous project:

- Files in ``checkpoints_parallel/*.done`` tell bcbio not to parallelize already
  executed pipeline tasks. This makes restarts faster by avoiding re-starting a
  cluster (when using distributed runs) for finished stages. If that behaviour
  is not desired for a task, removing the checkpoint file will get things
  parallelizing again.

- If the processing of a task is nearly finished the last jobs of this task will be
  running and bcbio will wait for those to finish.

IPython parallelization problems
================================

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
processes. The values are specified on a *memory-per-core* basis and
thus bcbio-nextgen handles memory scheduling by:

- :ref:`parallel-machine`

- Calculating the memory and core usage.
  The system configuration :ref:`config-resources` contains the
  expected core and memory usage of external programs.

- Adjusting the specified number of total cores to avoid
  over-scheduling memory. This allows running programs with more than
  the available memory per core without getting out of memory system
  errors.

- Passing total memory usage along to schedulers. The SLURM, SGE,
  Torque and PBSPro schedulers use this information to allocate memory to
  processes, avoiding issues with other scheduled programs using
  available memory on a shared machine.

As a result of these calculations, the cores used during processing
will not always correspond to the maximum cores provided in the input
``-n`` parameter. The goal is rather to intelligently maximize cores and
memory while staying within system resources. Note that memory
specifications are for a single core, and the pipeline takes care of
adjusting this to actual cores used during processing.

.. _parallel-machine:

Determining available cores and memory per machine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bcbio automatically tries to determine the total available memory and cores per
machine for balancing resource usage. For multicore runs, it retrieves total
memory from the current machine. For parallel runs, it spawns a job on the queue
and extracts the system information from that machine. This expects a
homogeneous set of machines within a cluster queue. You can see the determined
cores and total memory in ``provenance/system-ipython-queue.yaml``.

For heterogeneous clusters or other cases where bcbio does not correctly
identify available system resources, you can manually set the machine cores and
total memory in the ``resource`` section of either a sample YAML file
(:ref:`sample-resources`) or ``bcbio_system.yaml``::

    resources:
      machine:
        memory: 48.0
        cores: 16

The memory usage is total available on the machine in Gb, so this specifies that
individual machines have 48Gb of total memory and 16 cores.

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

In addition to open file handle limits (``ulimit -n``) large processes may also
run into issues with available max user processes (``ulimit -u``). Some systems
set a low soft limit (``ulimit -Su``) like 1024 but a higher hard limit
(``ulimit -Hu``), allowing adjustment without root privileges. The IPython
controllers and engines do this automatically, but the main ``bcbio_nextgen.py``
driver process cannot. If this scheduler puts this process on the same node as
worker processes, you may run into open file handle limits due to work happening
on the workers. To fix this, manually set ``ulimit -u a_high_number`` as part of
the submission process for the main process.

For a Ubuntu system, edit ``/etc/security/limits.conf`` to set the
soft and hard ``nofile`` descriptors, and edit
``/etc/pam.d/common-session`` to add ``pam_limits.so``. See
`this blog post`_ for more details.

For CentOS/RedHat systems, edit ``/etc/security/limits.conf`` and
``/etc/security/limits.d/90-nproc.conf`` to `increase maximum open files and
user limits <http://ithubinfo.blogspot.com/2013/07/how-to-increase-ulimit-open-file-and.html>`_.

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

- Harvard and Dell: See the 'Distributed File Systems' section of our
  `post on scaling bcbio-nextgen`_ for details about the setup within
  `Harvard FAS Research Computing`_ and thoughts on scaling and
  hardware. We also collaborate with Dell to
  test the pipeline on `Dell's Active Infrastructure for Life Sciences`_.
  We found the biggest initial factor limiting scaling was network
  bandwidth between compute and storage nodes.

.. _post on scaling bcbio-nextgen: http://bcb.io/2013/05/22/scaling-variant-detection-pipelines-for-whole-genome-sequencing-analysis/
.. _Harvard FAS Research Computing: http://rc.fas.harvard.edu/
.. _Dell's Active Infrastructure for Life Sciences: http://dell.com/ai-hpc-lifesciences

Spark
=====
Some GATK tools like recalibration use Apache Spark for parallelization. By default
bcbio runs these with multicore parallelization on a single node, to fit in standard
cluster and local compute environments. If you have a custom Spark cluster on your system
you can use that for GATK by setting up the appropriate configuration in your
:ref:`sample-resources`::

    resources:
        gatk-spark:
            options: [--spark-master, 'spark://your-spark-cluster:6311']
