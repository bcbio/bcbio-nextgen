.. _docs-cloud:

Amazon Web Services
-------------------

`Amazon Web Services (AWS) <https://aws.amazon.com/>`_ provides a flexible cloud
based environment for running analyses. Cloud approaches offer the ability to
perform analyses at scale with no investment in local hardware. They also offer
full programmatic control over the environment, allowing bcbio to automate the
entire setup, run and teardown process.

`bcbio-vm <https://github.com/chapmanb/bcbio-nextgen-vm>`_ provides a wrapper
around bcbio-nextgen that automates interaction with AWS and `Docker
<https://www.docker.com/>`_. ``bcbio_vm.py`` also cleans up the command line
usage to make it more intuitive and provides a superset of functionality
available in ``bcbio_nextgen.py``. bcbio-vm uses `Elasticluster
<https://github.com/gc3-uzh-ch/elasticluster>`_ to build a cluster on AWS with
an encrypted NFS mounted drive and an optional Lustre shared filesystem.

Local setup
===========

``bcbio_vm.py`` provides the automation to start up and administer remote bcbio
runs on AWS. This only requires a local installation of the python wrapper code,
not any of the Docker containers or biological data, which will all get
installed on AWS. The easier way to install is using `conda`_ with an isolated
Python::

    wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
    bash Miniconda-latest-Linux-x86_64.sh -b -p ~/install/bcbio-vm/anaconda
    ~/install/bcbio-vm/anaconda/bin/conda install --yes -c https://conda.binstar.org/bcbio bcbio-nextgen-vm
    ln -s ~/install/bcbio-vm/anaconda/bin/bcbio_vm.py /usr/local/bin/bcbio_vm.py

We support both Linux and Mac OSX as clients for running remote AWS bcbio clusters.

.. _conda: http://conda.pydata.org/

Data preparation
================

The easiest way to organize AWS projects is using an analysis folder inside an
`S3 bucket <http://aws.amazon.com/s3/>`_. Create a bucket and folder for your analysis and
upload fastq, BAM and, optionally, a region BED file. Bucket names should
include only lowercase letters, numbers and hyphens (``-``) to conform to
`S3 bucket naming restrictions <http://docs.aws.amazon.com/AmazonS3/latest/dev/BucketRestrictions.html>`_
and avoid issues with resolution of SSL keys. You can create buckets and upload
files using the
`AWS S3 web console <https://console.aws.amazon.com/s3/>`_,
`the AWS cli client <http://aws.amazon.com/cli/>`_ or specialized tools
like `gof3r <https://github.com/rlmcpherson/s3gof3r>`_.

You will also need a template file describing the type of run to do and a CSV
file mapping samples in the bucket to names and any other metadata. See the
:ref:`automated-sample-config` docs for more details about these files. Also
upload both of these files to S3.

With that in place, prepare and upload the final configuration to S3 with::

    bcbio_vm.py template s3://your-project/your-analysis/template.yaml s3://your-project/your-analysis/name.csv

This will find the input files in the ``s3://your-project/your-analysis`` bucket, associate
fastq and BAM files with the right samples, and add a found BED files as
``variant_regions`` in the configuration. It will then upload the final
configuration back to S3 as ``s3://your-project/your-analysis/name.yaml``, which you can run
directly from a bcbio cluster on AWS.

We currently support human analysis with both the GRCh37 and hg19 genomes. We
can also add additional genomes as needed by the community and generally welcome
feedback and comments on reference data support.

Extra software
~~~~~~~~~~~~~~

We're not able to automatically install some useful tools in pre-built docker
containers due to licensing restrictions. Variant calling with GATK requires a
manual download from the `GATK download`_ site for academic users.  Appistry
provides `a distribution of GATK for commercial users`_. Commercial users also
need a license for somatic calling with muTect. To make these jars available,
upload them to the S3 bucket in a ``jars`` directory. bcbio will automatically
include the correct GATK and muTect directives during your run.  Alternatively,
you can also manually specify the path to the jars using the global
``resources`` section of your input sample YAML file::

    resources:
      gatk:
        jar: s3://bcbio-syn3-eval/jars/GenomeAnalysisTK.jar

.. _GATK download: http://www.broadinstitute.org/gatk/download
.. _a distribution of GATK for commercial users: http://www.appistry.com/gatk

AWS setup
=========

The first time running bcbio on AWS you'll need to setup permissions, VPCs and
local configuration files. We provide commands to automate all these steps and once
finished, they can be re-used for subsequent runs. To start you'll need to have
an account at Amazon and your Access Key ID and Secret Key ID from the
`AWS security credentials page
<https://console.aws.amazon.com/iam/home?#security_credential>`_. These can be
`IAM credentials <https://aws.amazon.com/iam/getting-started/>`_ instead of root
credentials as long as they have administrator privileges. Make them available
to bcbio using the standard environmental variables::

  export AWS_ACCESS_KEY_ID=your_access_key
  export AWS_SECRET_ACCESS_KEY=your_secret_key

With this in place, two commands setup your elasticluster and AWS environment to
run a bcbio cluster. The first creates public/private keys, a bcbio IAM user,
and sets up an elasticluster config in ``~/.bcbio/elasticluster/config``::

  bcbio_vm.py aws iam

The second configures a VPC to host bcbio::

  bcbio_vm.py aws vpc

The ``aws vpc`` command is idempotent and can run multiple times if you change or
remove parts of the infrastructure. You can also rerun the ``aws iam`` command,
but if you'd like to generate a new elasticluster configuration file
(``~/.bcbio/elasticluster/config``) add the recreate flag: ``bcbio_vm.py aws iam
--recreate``. This generates a new set of IAM credentials and public/private
keys. These are only stored in the ``~/.bcbio`` directory so you need to fully
recreate them if you delete the old ones.

Running a cluster
=================

Following this setup, you're ready to run a bcbio cluster on AWS. We start
from a standard Ubuntu AMI, installing all software for bcbio and the cluster as
part of the boot process.

To configure your cluster run::

   bcbio_vm.py aws config edit

This dialog allows you to define the cluster size and machine resources you'd
like to use. The defaults only have small instances to prevent accidentally
starting an `expensive run <http://aws.amazon.com/ec2/pricing/>`_. If you're
planning a run with less than 32 cores, do not use a cluster and instead run
directly on a single machine using one of the `large r3 or c3 instances
<http://aws.amazon.com/ec2/instance-types/>`_.

This script also sets the size of the `encrypted NFS-mounted drive
<http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/EBSEncryption.html>`_, which
you can use to store processing data when running across a distributed
cluster. At scale, you can replace this with a Lustre shared filesystem. See
below for details on launching and attaching a Lustre filesystem to a cluster.

To ensure everything is correctly configured, run::

    bcbio_vm.py aws info

When happy with your setup, start the cluster with::

    bcbio_vm.py aws cluster start

The cluster will take five to ten minutes to start. If you encounter any
intermittent failures, you can rerun the cluster configuration step with
``bcbio_vm.py aws cluster setup`` or the bcbio-specific installation with
``bcbio_vm.py aws cluster bootstrap``.

Running Lustre
==============

Elasticluster mounts the ``/encrypted`` directory as a NFS share available
across all of the worker machines. You can use this as a processing directory
for smaller runs but for larger runs may need a scalable distributed file
system. bcbio supports using
`Intel Cloud Edition for Lustre (ICEL) <https://wiki.hpdd.intel.com/display/PUB/Intel+Cloud+Edition+for+Lustre*+Software>`_
to set up a Lustre scratch filesystem on AWS.

- Subscribe to `ICEL in the Amazon Marketplace
  <https://aws.amazon.com/marketplace/pp/B00GK6D19A>`_.

- By default, the Lustre filesystem will be 2TB and will be accessible to
  all hosts in the VPC. Creation takes about ten minutes and can happen in
  parallel while elasticluster sets up the cluster. Start the stack::

    bcbio_vm.py aws icel create

  If you encounter any intermittent failures when installing collectl plugin, that
  means lustre server is created successfully, you can rerun the lustre configuration step
  with ``bcbio_vm.py aws icel create --setup``. If you had any failure creating the lustre
  server before the collectl plugin installation, you should stop it, and try again.
  

- Once the ICEL stack and elasticluster cluster are both running, mount the
  filesystem on the cluster::

    bcbio_vm.py aws icel mount

- The cluster instances will reboot with the Lustre filesystem mounted.

Running an analysis
===================

To run the analysis, connect to the head node with::

    bcbio_vm.py aws cluster ssh

If you started a single machine, or a cluster using encrypted NFS, run with::

    mkdir /encrypted/your-project
    cd !$ && mkdir work && cd work
    bcbio_vm.py run -n 8 s3://your-project/your-analysis/name.yaml

Where the ``-n`` argument should be the number of cores on the machine.

To run on a full cluster with a Lustre filesystem::

    sudo mkdir /scratch/cancer-dream-syn3-exome
    sudo chown ubuntu !$
    cd !$ && mkdir work && cd work
    bcbio_vm.py ipythonprep s3://your-project/your-analysis/name.yaml slurm cloud -n 60
    sbatch bcbio_submit.sh

Where 60 is the total number of cores to use across all the worker nodes.  Of
your total machine cores, allocate 2 for the base bcbio_vm script and IPython
controller instances. The `SLURM workload manager <http://slurm.schedmd.com/>`_
distributes jobs across your cluster on a queue called ``cloud``.  A
``slurm-PID.out`` file in the work directory contains the current status of the
job, and ``sacct_std`` provides the status of jobs on the cluster. If you are
new to SLURM, here is a summary of useful
`SLURM commands <https://rc.fas.harvard.edu/resources/running-jobs/#Summary_of_SLURM_commands>`_.

On successful completion, bcbio uploads the results of the analysis back into your s3
bucket and folder as ``s3://your-project/your-analysis/final``. You can now cleanup the cluster and
Lustre filesystem.

Graphing resource usage
=======================

AWS runs include automatic monitoring of resource usage with
`collectl <http://collectl.sourceforge.net/>`_. bcbio_vm uses collectl statistics
to plot CPU, memory, disk and network usage during each step of a run. To
prepare resource usage plots after finishing an analysis, first copy the
``bcbio-nextgen.log`` file to your local computer. Either use
``bcbio_vm.py elasticluster sftp bcbio`` to copy from the work directory on AWS
(``/encrypted/your-project/work/log/bcbio-nextgen.log``) or transfer it from the
output S3 bucket (``your-project/your-analysis/final/DATE_your-project/bcbio-nextgen.log``).

If your run worked cleanly you can use the log input file directly. If you had
failures and restarts, or would only like to graph part of the run, you can edit
the timing steps. Run ``grep Timing bcbio-nextgen.log > your-run.txt`` to get
the timing steps only, then edit as desired.

Retrieve the collectl statistics from the AWS cluster and prepare the resource
usage graphs with::

    bcbio_vm.py graph bcbio-nextgen.log

Collectl stats will be in ``monitoring/collectl`` and plots are in
``monitoring/graphs``. If you need to re-run plots later after shutting the
cluster down, you can use the local collectl stats instead of retrieving from
the server by running ``bcbio_vm.py graph bcbio-nextgen.log --cluster none``.
If you'd like to run graphing from a local non-AWS run, manually place collectl
files from each node to analyze in ``monitoring/collectl/yournodename-timestamp.raw.gz``.
In addition to plots, the
`summarize_timing.py <https://github.com/chapmanb/bcbio-nextgen/blob/master/scripts/utils/summarize_timing.py>`_
utility script prepares a summary table of run times per step.

Shutting down
=============

The bcbio Elasticluster and Lustre integration can spin up a lot of AWS
resources. You'll be paying for these by the hour so you want to clean them up
when you finish running your analysis. To stop the cluster::

    bcbio_vm.py aws cluster stop

To remove the Lustre stack::

    bcbio_vm.py aws icel stop

Double check that all instances have been properly stopped by looking in the AWS
console.

Manual configuration
====================

Experienced `elasticluster <https://github.com/gc3-uzh-ch/elasticluster>`_ users
can edit the configuration files themselves. bcbio provides a small wrapper
that automatically reads and writes these configurations to avoid users needing
to understand elasticluster internals, but all functionality is fully available.
Edit your ``~/.bcbio/elasticluster/config`` file to change parameters. You can
also see the `latest example configuration
<https://github.com/chapmanb/bcbio-nextgen-vm/blob/master/elasticluster/config>`_.
in the bcbio-vm GitHub repository for more details on the other available options.
