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

    bcbio_nextgen.py bcbio_system.yaml bcbio_sample.yaml -t local -n 12

IPython parallel
~~~~~~~~~~~~~~~~

`IPython parallel`_ provides a distributed framework for performing
parallel computation in standard cluster environments. The
bcbio-nextgen setup script installs both IPython and `pyzmq`_, which
provides Python bindings for the `ZeroMQ`_ messaging library.

Run an analysis using ipython for parallel execution::

    bcbio_nextgen.py bcbio_system.yaml bcbio_sample.yaml -t ipython -n 12 -s lsf -q queue

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

.. _IPython parallel: http://ipython.org/ipython-doc/dev/index.html
.. _pyzmq: https://github.com/zeromq/pyzmq
.. _ZeroMQ: http://www.zeromq.org/
