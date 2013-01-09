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

We also provide a setup script that correctly configures IPython for
different cluster environments. Run this once on new machines::

    bcbio_nextgen_setup.py -s lsf -q your_queue_name

The ``-s`` flag specifies a type of scheduler to use ``(lsf, sge)``.
The ``-q`` flag specifies the queue to submit jobs to. The
setup script will create IPython configurations for parallel and
multicore jobs.

When setup, run an analysis specifying ipython for parallel execution::

    bcbio_nextgen.py bcbio_system.yaml bcbio_sample.yaml -t ipython -n 12

Celery and RabbitMQ
~~~~~~~~~~~~~~~~~~~

We still support celery and RabbitMQ messaging, but please try IPython
when setting up a new cluster. The IPython approach is under active
development and supports additional cluster and parallel approaches.

To enable parallel messaging:

1. Configure RabbitMQ as described below. Ensure all processing machines
   can talk to the RabbitMQ server on port 5672. Update
   ``universe_wsgi.ini`` to contain the server details.

2. Edit your ``post_process.yaml`` file to set parameters in the
   ``distributed`` section corresponding to your environment: this
   includes the type of cluster management and arguments to start jobs.

3. Run ``bcbio_nextgen.py`` with parameters for a distributed cluster
   environment. It takes care of starting worker nodes, running the
   processing, and then cleaning up after jobs::

      bcbio_nextgen.py post_process.yaml flowcell_dir run_info.yaml
                       -t messaging -n 20

RabbitMQ configuration
**********************

RabbitMQ messaging manages communication between the sequencing machine
and the analysis machine. This allows complete separation between all of
the machines. The RabbitMQ server can run anywhere; an easy solution is
to install it on the Galaxy and analysis server::

        (yum or apt-get) install rabbitmq-server

Setup rabbitmq for passing Galaxy and processing messages::

        rabbitmqctl add_user <username> <password>
        rabbitmqctl add_vhost bionextgen
        rabbitmqctl set_permissions -p bionextgen <username> ".*" ".*" ".*"

Then adjust the ``[galaxy_amqp]`` section of your ``universe_wsgi.ini``
Galaxy configuration file. An example configuration is available in the
config directory; you'll need to specifically change these three values::

        [galaxy_amqp]
        host = <host you installed the RabbitMQ server on>
        userid = <username>
        password = <password>

ssh keys
********

The sequencing, analysis and storage machines transfer files using
secure copy. This requires that you can securely copy files between
machines without passwords, using `ssh public key`_ authentication.
You want to enable password-less ssh for the following machine
combinations:

-  Analysis server to ``illumina_finished_msg`` machine
-  Storage server to ``illumina_finished_msg`` machine

Sequencing machines
*******************

The sequencer automation has been fully tested using Illumina GAII and
HiSeq sequencing machines. The framework is general and supports other
platforms; we welcome feedback from researchers with different machines
at their institutions.

Illumina machines produce run directories that include the date, machine
identifier, and flowcell ID::

    110216_HWI-EAS264_00035_FC638NPAAXX

A shortened name, with just date and flowcell ID, is used to uniquely
identify each flowcell during processing.

.. _ssh public key: http://macnugget.org/projects/publickeys/
.. _IPython parallel: http://ipython.org/ipython-doc/dev/index.html
.. _pyzmq: https://github.com/zeromq/pyzmq
.. _ZeroMQ: http://www.zeromq.org/
