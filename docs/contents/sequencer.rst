Sequencer integration
---------------------

bcbio-nextgen supports processing of samples arriving from sequencing
machines. It automates the generation of de-multiplexed fastq files and
sample configurations that feed directly into standard bcbio-nextgen
pipelines.

This requires an Illumina sequencer, and samples annotated using the basic
`Galaxy nglims`_. The approach is general and we'd be happy to
collaborate or accept contributions for supporting other sequencers or LIMS
systems.

.. _Galaxy nglims: https://wiki.galaxyproject.org/Admin/SampleTracking/NextGen

Overview
********

A fully automated setup consists of three components:

- A front end `Galaxy nglims`_ system. Users provide details on samples,
  including organism information and choice of post-processing
  pipeline. Sequencing operators annotate multiplexed flowcells with the
  locations of samples. The automated processing scripts use this information to
  prepare sample configurations for downstream processing.

- A sequencer output machine, where the sequencer dumps reads. A
  cronjob detects newly finished flowcells and processes output into
  demultiplexed fastq.

- An analysis machine where bcbio-nextgen analysis occurs. The cronjob on the
  sequencer output machines transfers fastq files and initiates multi-core
  processing. On completion, the analysis machine uploads results to an attached
  Galaxy instance.

Sequencer output machine
************************

Post-sequencing processing, including demultiplexing, initiate via a cronjob run
on the Illumina output machine::

    PATH=/usr/local/bin
    @hourly bcbio_nextgen.py sequencer /opt/bcbio/transfer_info.yaml

`transfer_info.yaml`_ is a configuration file specifying locations of output
directories where sequencer results appear. It also contains the location of the
Galaxy nglims server to retrieve sample details, and the downstream analysis
server to transfer files and start bcbio-nextgen pipelines.

Illumina machines produce run directories that include the date, machine
identifier, and flowcell ID::

    110216_HWI-EAS264_00035_FC638NPAAXX

the cronjob script identifies directories with newly finished
output. Unprocessed directories, identified by the date and flowcell ID,
continue on for demultiplexing, transfer and analysis.

Illumina's `bcl2fastq`_ performs demultiplexing and conversion to fastq. The
``configureBclToFastq.pl`` script from `bcl2fastq`_ needs to be available on the PATH
specified within your cronjob.

Following demultiplexing, the script combines separate fastq files into single
or paired fastq files per sample. This includes combining samples multiplexed
across multiple lanes in a flowcell; any samples with the same
project and sample name get combined.

Finally the script prepares a sample configuration file defining processing
steps to perform, transfers the configuration and fastq files to the analysis
machine, and initiates processing. Transfer between the sequencer output and
analysis machines occurs using secure rsync, which requires the ability to
securely login between machines without passwords using `ssh public key`_
authentication.

.. _transfer_info.yaml: https://github.com/chapmanb/bcbio-nextgen/blob/master/config/transfer_info.yaml
.. _bcl2fastq: http://support.illumina.com/downloads/bcl2fastq_conversion_software_184.ilmn
.. _ssh public key: http://macnugget.org/projects/publickeys/

Analysis server
***************

The analysis server runs bcbio-nextgen pipelines and uploads results to a local
Galaxy server. A bcbio-nextgen server receives processing commands and start
them. The following command starts a server on port 8080 using 16 cores for analysis::

     bcbio_nextgen.py server p 8080 -n 16

To make this available outside of the current machine use a proxy server like
`nginx`_ with the following configuration::

    upstream bcbio {
      server localhost:8080;
    }

    server {
        listen       80 default_server;
        location /bcbio/ {
           proxy_pass http://bcbio/;
           proxy_set_header   X-Forwarded-Host $host;
           proxy_set_header   X-Forwarded-For  $proxy_add_x_forwarded_for;
        }
    }

.. _nginx: http://nginx.org/

Specify this URL as ``server: http://analysis.you.org/bcbio`` in your
`transfer_info.yaml`_ file to enable the sequencing output machine to
communicate with the analysis server.

The analysis server currently handles multicore processing. We'd be happy to
collaborate on approaches to allow it to automatically start bcbio-nextgen jobs
on HPC clusters or other types of distributed environments.

On analysis completion the pipeline transfers processed files to a Galaxy server
based on a pre-specified :ref:`upload-configuration` configuration. This
includes the alignment BAM files, quality control, and other pipeline
specific files like variant calls or RNA-seq counts. It organizes files in
project and sample specific folders within Galaxy's data libraries, making them
available to researchers for additional analysis.

Debugging
*********

This section contains tips and tricks on restarting processing in case of
problems. Flowcell processing occurs under the directory specified by
``process -> dir`` in your `transfer_info.yaml`_ file. Each flowcell directory
contains the sample YAML configuration file, an analysis directory,
and demultiplexed fastqs::

    140313_SN728_0206_AC3KL2ACXX
    ├── 140313_SN728_0206_AC3KL2ACXX.csv
    ├── 140313_SN728_0206_AC3KL2ACXX.yaml
    ├── analysis
    ├── Data
    ├── fastq
    ├── RunInfo.xml
    └── runParameters.xml

The full log file for a processing run is in
``analysis/log/bcbio-nextgen-debug.log`` and will contain useful details about
why a run failed. You can manually restart processing of a run using the
standard bcbio-nextgen command line::

    cd FLOWCELL/analysis
    bcbio_nextgen.py ../../FLOWCELL ../FLOWCELL.yaml -n 16
