## Overview

Python scripts and modules for automated next gen sequencing analysis.
These provide a fully automated pipeline for taking sequencing results
from an Illumina sequencer, converting them to standard Fastq format,
aligning to a reference genome, doing SNP calling, and producing a
summary PDF of results.

The scripts are tightly integrated with the [Galaxy][1] 
web-based analysis tool. Samples are entered and tracked through a LIMS
system and processed results are uploading into Galaxy Data Libraries for
researcher access and additional analysis. Our [clone of Galaxy][2]
tracks the main development work adding an intuitive interface for sample
management on top of the existing functionality.

[1]: http://galaxy.psu.edu/
[2]: http://bitbucket.org/chapmanb/galaxy-central

ToDo: pictoric representation of the workflow/pipeline

## Code structure

The main scripts that handle automation of the analysis and storage
are:

* `scripts/illumina_finished_msg.py` -- Sits on a machine where sequencing
  runs are dumped. Should be run via a cron job every hour or so to
  check for new output runs.

* `scripts/analyze_finished_sqn.py` -- Continuously running server script on
  the Galaxy analysis machine. When new results are reported in the messaging queue,
  this copies over the relevant fastq files and starts an automated
  analysis.

* `scripts/store_finished_sqn.py` -- Continuously running server that
  manages long term storage of larger files, like qseq and images.

System specific information is specified in configuration files:

* `config/transfer_info.yaml` -- Configuration on the sequencing
  machine, specifying where to check for new sequencing data.
* `config/post_process.yaml` -- Configuration for analysis and
  storage. This contains links to Galaxy, program commandlines and
  customization for processing algorithms.
* `config/universe_wsgi.ini` -- Configuration variables that should be
  included on your Galaxy server. RabbitMQ details are specified here
  for communication between the sequencing and analysis machine.

The scripts involved in the actual processing:

* `scripts/automated_initial_analysis.py` -- Drives the high level analysis of
  sequencing lanes based on information specified through the Galaxy LIMS system
* `scripts/upload_to_galaxy.py` -- Handles storing and uploading Fastq,
  alignment, analysis and summary files to Galaxy.
* `scripts/align_summary_report.py` -- Produces a PDF summary file with
  statistics on alignments, duplicates, GC distribution, quality scores,
  and other metrics of interest.

## Installation

### Backend

Clone a copy from chapmanb branch:

git clone git://github.com/chapmanb/bcbb.git

Install using the standard python approach. You will need a recent
version of Python 2 (2.6 or 2.7). The required python library
dependencies will also be installed:

        cd bcbb/nextgen && python setup.py install

Copy the YAML & ini files in config and adjust them to match your
environment. It is also a good idea to set your $PATH pointing to
any third-party binaries you are using.

#### RabbitMQ messaging server

Communication between the sequencing machine and the analysis machine
running Galaxy is managed using RabbitMQ messaging. This allows
complete separation between all of the machines. The RabbitMQ server
can run anywhere; an easy solution is to install it on the Galaxy
analysis server:

        (yum or apt-get) install rabbitmq-server

Setup rabbitmq for passing Galaxy messages:

        rabbitmqctl add_user <username> <password>
        rabbitmqctl add_vhost galaxy_messaging_engine
        rabbitmqctl set_permissions -p galaxy_messaging_engine <username> ".*" ".*" ".*"

Then adjust the `[galaxy_amqp]` section of your `universe_wsgi.ini`
Galaxy configuration file. An example configuration is available in
the config directory; you'll need to specifically change these three
values:

        [galaxy_amqp]
        host = <host you installed the RabbitMQ server on>
        userid = <username>
        password = <password>

#### ssh keys

The sequencing, analysis and storage machines transfer files using
secure copy. This requires that you can securely copy files between
machines without passwords, using [ssh public key authentication][i1].
You want to enable password-less ssh for the following machine
combinations:

* `analyze_finished_sqn` server to `illumina_finished_msg` machine
* `store_finished_sqn` server to `illumina_finished_msg` machine
* `analyze_finished_sqn` to actual analysis machine, if different else
  localhost
* `store_finished_sqn` server to actual storage machine, if different
  else localhost

[i1]: http://macnugget.org/projects/publickeys/

### Frontend

Follow up the following instructions to setup the Galaxy instance:

https://bitbucket.org/galaxy/galaxy-central/wiki/LIMS/nglims

### Development environment

The installation instructions assume that you have full root access to install
python modules and packages (production environment). If this is not the case,
you may want to install a python VirtualEnv and other tools automatically on your $HOME
to ease your development needs using the following script:

wget https://bitbucket.org/brainstorm/custom_env/raw/1cd4f4ae27d5/pyHost.sh && ./pyHost.sh

## Requirements

### Next gen analysis

* [Picard][3]
* [GATK][4]
* [bowtie][5]
* [fastx toolkit][6]
* [SolexaQA][6b]
* [samtools][7]
* [snpEff][16]

[3]: http://picard.sourceforge.net/
[4]: http://www.broadinstitute.org/gsa/wiki/index.php/The_Genome_Analysis_Toolkit
[5]: http://bowtie-bio.sourceforge.net/
[6]: http://hannonlab.cshl.edu/fastx_toolkit/
[6b]: http://solexaqa.sourceforge.net/
[7]: http://samtools.sourceforge.net/
[16]: http://sourceforge.net/projects/snpeff/

### Processing infrastructure

* RabbitMQ
* LaTeX -- pdflatex
* R with ggplot2 and sqldf
* ps2pdf

### Python modules

* [Biopython][10]
* [rpy2][11]
* [pysam][12]
* [mako][13]
* [PyYAML][14]
* [amqplib][15]
* [logbook] [17]

[10]: http://biopython.org
[11]: http://rpy.sourceforge.net/rpy2.html
[12]: http://code.google.com/p/pysam/
[13]: http://www.makotemplates.org/
[14]: http://pyyaml.org/
[15]: http://code.google.com/p/py-amqplib
[17]: http://packages.python.org/Logbook
