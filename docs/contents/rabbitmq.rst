Messaging
---------

RabbitMQ
~~~~~~~~~

RabbitMQ messaging manages communication between the sequencing machine
and the analysis machine. This allows complete separation between all of
the machines. The RabbitMQ server can run anywhere; an easy solution is
to install it on the Galaxy and analysis server:

::

        (yum or apt-get) install rabbitmq-server

Setup rabbitmq for passing Galaxy and processing messages:

::

        rabbitmqctl add_user <username> <password>
        rabbitmqctl add_vhost bionextgen
        rabbitmqctl set_permissions -p bionextgen <username> ".*" ".*" ".*"

Then adjust the ``[galaxy_amqp]`` section of your ``universe_wsgi.ini``
Galaxy configuration file. An example configuration is available in the
config directory; you'll need to specifically change these three values:

::

        [galaxy_amqp]
        host = <host you installed the RabbitMQ server on>
        userid = <username>
        password = <password>

ssh keys
~~~~~~~~

The sequencing, analysis and storage machines transfer files using
secure copy. This requires that you can securely copy files between
machines without passwords, using `ssh public key`_ authentication.
You want to enable password-less ssh for the following machine
combinations:

-  Analysis server to ``illumina_finished_msg`` machine
-  Storage server to ``illumina_finished_msg`` machine

Sequencing machines
~~~~~~~~~~~~~~~~~~~

The sequencer automation has been fully tested using Illumina GAII and
HiSeq sequencing machines. The framework is general and supports other
platforms; we welcome feedback from researchers with different machines
at their institutions.

Illumina machines produce run directories that include the date, machine
identifier, and flowcell ID:

::

    110216_HWI-EAS264_00035_FC638NPAAXX

A shortened name, with just date and flowcell ID, is used to uniquely
identify each flowcell during processing.

.. _ssh public key: http://macnugget.org/projects/publickeys/
