Deprecated functionality
------------------------

This section describes older functionality migrated to new approaches.
We maintain support for back-compatibility purposes but suggest moving
to the updated in-development approaches.

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

Sequencer support
~~~~~~~~~~~~~~~~~

The main scripts that handle automation of the analysis and storage are:

-  ``scripts/illumina_finished_msg.py`` -- Sits on the sequencer output
   machine; run via a cron job every hour to check for new output runs.

-  ``scripts/nextgen_analysis_server.py`` -- Main server script; runs
   specific servers for top level analysis management, storage, or
   distributed processing of individual analysis steps:

-  Called with ``-q toplevel`` -- Coordinate sample processing, running
   the full automated analysis and optionally uploading results to
   Galaxy. ``illumina_finished_msg.py`` and a Galaxy graphical front end
   both send messages on this queue for starting processing.

-  Called with no queue arguments -- Run individual steps in the
   analysis pipeline. Start multiple servers on distributed machines
   connected with a shared filesystem to allow scaling on a cluster or
   Amazon EC2.

-  Called with ``-q storage`` -- Manage long term storage of larger
   files, like qseq and images.

Specify system specific information in configuration files:

-  ``config/transfer_info.yaml`` -- Configuration on the sequencing
   machine, specifying where to check for new sequencing data.
-  ``config/post_process.yaml`` -- Configuration for analysis and
   storage. This contains links to Galaxy, program locations and
   customization for processing algorithms.
-  ``config/universe_wsgi.ini`` -- Variables used from your Galaxy
   server configuration, including RabbitMQ details for communication
   between the sequencing and analysis machines.
