.. _docs-cloud:

Cloud
-----

bcbio has two approaches to running on cloud providers like
`Amazon Web Services (AWS) <https://aws.amazon.com/>`_,
`Google Cloud (GCP) <https://cloud.google.com/>`_ and
`Microsoft Azure <https://azure.microsoft.com>`_. For smaller projects
we use a `simplified ansible based
approach
<https://github.com/bcbio/bcbio-nextgen/tree/master/scripts/ansible#simplified-bcbio-cloud-usage>`_
which automates spinning up single multicore machines for running either
traditional or :ref:`docs-cwl` bcbio runs.

For larger distributed projects, we're actively working on using :ref:`docs-cwl`
support with runners like `Cromwell <http://cromwell.readthedocs.io>`_ that
directly interface and run on cloud services. We'll document these approaches
here as they're tested and available.

For getting started, the CWL :ref:`docs-cwl-installation` documentation
describes how to install `bcbio-vm <https://github.com/bcbio/bcbio-nextgen-vm>`_,
which provides a wrapper around bcbio that automates interaction with cloud
providers and `Docker <https://www.docker.com/>`_. ``bcbio_vm.py`` also cleans
up the command line usage to make it more intuitive and provides a superset of
functionality available in ``bcbio_nextgen.py``.

.. _docs-cloud-gcp:

Google Cloud (GCP)
##################

Cromwell runs bcbio CWL pipelines on Google Cloud using the
`Google Pipelines API <https://cloud.google.com/genomics/reference/rest/>`_.

Setup
=====

To setup a Google Compute environment, you'll make use of the `Web based console
<https://console.cloud.google.com>`_ and `gcloud and gsutil from the Google
Cloud SDK <https://cloud.google.com/sdk/>`_, which provide command line
interfacts to manage data in Google Storage and Google Compute instances. You
can install with::

    bcbio_conda install -c conda-forge -c bioconda google-cloud-sdk

For authentication, you want to set up a `Google Cloud Platform service account
<https://cloud.google.com/docs/authentication/production>`_. The environmental variable
``GOOGLE_APPLICATION_CREDENTIALS`` identifies a
`JSON file of credentials <https://cloud.google.com/docs/authentication/getting-started>`_.
which bcbio passes to Cromwell for authentication::

    gcloud auth login
    gcloud projects create your-project
    gcloud iam service-accounts create your-service-account
    gcloud projects add-iam-policy-binding your-project --member \
      "serviceAccount:your-service-account@your-project.iam.gserviceaccount.com" --role "roles/owner"
    gcloud iam service-accounts keys create ~/.config/gcloud/your-service-account.json \
      --iam-account your-service-account@your-project.iam.gserviceaccount.com
    export GOOGLE_APPLICATION_CREDENTIALS=~/.config/gcloud/your-service-account.json

You'll need a project for your run along, with the Google Genomics API enabled,
and a Google Storage bucket for your data and run intermediates::

    gcloud config set project your-project
    gcloud services enable genomics.googleapis.com
    gsutil mb gs://your-project

Additional documentation for Cromwell: `Google Pipelines API
<https://cromwell.readthedocs.io/en/stable/tutorials/PipelinesApi101/>`_ and
`Google authentication <https://github.com/broadinstitute/cromwell/blob/develop/docs/backends/Google.md>`_.

Data preparation
================

Cromwell can localize data present in Google Storage buckets as part of the run
process and bcbio will translate the data present in these storage bucket into
references for the CWL run inputs.

Upload your data with ``gsutil``::

    gsutil cp your_data.bam gs://your-project/inputs/


Create a ``bcbio_system-gcp.yaml`` input file for :ref:`docs-cwl-generate`::

    gs:
      ref: gs://bcbiodata/collections
      inputs:
        - gs://your-project/inputs
    resources:
      default: {cores: 8, memory: 3G, jvm_opts: [-Xms750m, -Xmx3000m]}

Then create a sample input CSV and template YAML file for
:ref:`automated-sample-config`. The first column of the CSV file should contain
references to your input files (``your_file.bam`` or
``your_file_R1.fastq.gz;your_file_R2.fastq.gz``), which avoids needing to specify the
inputs on the command line.

Generate a Common Workflow Language representation::

   bcbio_vm.py template --systemconfig bcbio_system-gcp.yaml ${TEMPLATE}-template.yaml $PNAME.csv
   bcbio_vm.py cwl --systemconfig bcbio_system-gcp.yaml $PNAME/config/$PNAME.yaml

Running
=======

Run the CWL using Cromwell by specifying the project and root Google Storage
bucket for intermediates::

    bcbio_vm.py cwlrun cromwell $PNAME-workflow --cloud-project your-project \
        --cloud-root gs://your-project/work_cromwell

Amazon Web Services (AWS Batch)
###############################

We're working to support `Amazon Web Services (AWS) <https://aws.amazon.com/>`_
using AWS Batch and Cromwell, following the `AWS for Genomics documentation
<https://docs.opendata.aws/genomics-workflows/>`_. This documents the current
work in progress; it is not yet fully running and needs
`additional Cromwell development <https://github.com/broadinstitute/cromwell/issues/4586>`_
for AWS CWL support.

Setup
=====

0. Optionally, create a bcbio `IAM user <https://aws.amazon.com/iam/>`_ and
   bcbio keypair for creating AWS Batch specific resources. bcbio-vm can
   automate this process, although they can also be pre-existing. If you'd like
   to use bcbio-vm automation, you'll need to have
   an account at Amazon and your Access Key ID and Secret Key ID from the
   `AWS security credentials page
   <https://console.aws.amazon.com/iam/home?#security_credential>`_. These can be
   `IAM credentials <https://aws.amazon.com/iam/getting-started/>`_ instead of root
   credentials as long as they have administrator privileges. Make them available
   to bcbio using the standard environmental variables::

       export AWS_ACCESS_KEY_ID=your_access_key
       export AWS_SECRET_ACCESS_KEY=your_secret_key

   With this in place, ceate public/private keys and a bcbio IAM user with::

       bcbio_vm.py aws iam --region=us-east-1

1. Use either existing credentials or those created by bcbio, setup `AWS Credentials
   <https://boto3.amazonaws.com/v1/documentation/api/latest/guide/quickstart.html#configuration>`_
   for accessing AWS resources from your machine by editing `~/.aws/credentials`::

       [default]
       aws_access_key_id = YOURACCESSID
       aws_secret_access_key = yoursecretkey
       region = us-east-1

2. Automation creation of resources for AWS Batch. This includes creating
   a `custom Amazon Machine Image (AMI) for AWS Batch
   <https://docs.opendata.aws/genomics-workflows/aws-batch/create-custom-ami/>`_,
   which allows automatic allocation of additional disk space during workflow
   runs. It also sets up an `AWS Batch environment, VPC and IAM for running workflows
   <https://docs.opendata.aws/genomics-workflows/aws-batch/configure-aws-batch-cfn/>`_.
   A single bcbio-vm commands runs both CloudFormation scripts::

       bcbio_vm.py aws cromwell --keypair bcbio --bucket bcbio-batch-cromwell-test

   This will output the S3 bucket and job queue for running Cromwell::

      AMI: ami-00bd75374ccaa1fc6
      Region: us-east-1
      S3 bucket: s3://your-project
      Job Queue (Spot instances): arn:aws:batch:us-east-1:678711657553:job-queue/GenomicsDefaultQueue-358a1deb9f4536b
      High priority Job Queue: arn:aws:batch:us-east-1:678711657553:job-queue/GenomicsHighPriorityQue-3bff21e3c4f44d4

Data preparation
================

The easiest way to organize AWS projects is using an analysis folder inside an
`S3 bucket <http://aws.amazon.com/s3/>`_. Create a bucket and folder for your analysis and
upload input files (fastq or BAM) and other associated files.. Bucket names should
include only lowercase letters, numbers and hyphens (``-``) to conform to
`S3 bucket naming restrictions <http://docs.aws.amazon.com/AmazonS3/latest/dev/BucketRestrictions.html>`_
and avoid issues with resolution of SSL keys. You can create buckets and upload
files using the `the AWS cli client <http://aws.amazon.com/cli/>`_ or
`AWS S3 web console <https://console.aws.amazon.com/s3/>`_::

    aws s3 sync /local/inputs s3://your-bucket/inputs

Create a ``bcbio_system-aws.yaml`` input file for :ref:`docs-cwl-generate`::

    s3:
      ref: s3://bcbiodata/collections
      inputs:
        - s3://your-bucket/inputs
    resources:
      default: {cores: 8, memory: 3G, jvm_opts: [-Xms750m, -Xmx3000m]}

Generate a Common Workflow Language representation::

   CLOUD=aws
   bcbio_vm.py template --systemconfig bcbio_system-$CLOUD.yaml ${TEMPLATE}-template.yaml $PNAME.csv
   bcbio_vm.py cwl --systemconfig bcbio_system-$CLOUD.yaml $PNAME/config/$PNAME.yaml

Running
=======

Run the CWL using Cromwell by specifying the batch job queue
`Amazon Resource Name (ARN) <https://docs.aws.amazon.com/general/latest/gr/aws-arns-and-namespaces.html>`_
and bucket from the setup process::

    bcbio_vm.py cwlrun cromwell $PNAME-workflow \
      -cloud-project arn:aws:batch:us-east-1:678711657553:job-queue/GenomicsDefaultQueue-358a1deb9f4536b \
      -cloud-root s3://your-project

Amazon Web Services (old)
#########################

We're phasing out this approach to AWS support in bcbio and are actively
moving to Common Workflow Language based approaches. This documents the old
`Elasticluster
<https://github.com/gc3-uzh-ch/elasticluster>`_ approach to build a cluster on AWS with
an encrypted NFS mounted drive and an optional Lustre shared filesystem.

Data preparation
================

You need a template file describing the type of run to do and a CSV
file mapping samples in the bucket to names and any other metadata. See the
:ref:`automated-sample-config` docs for more details about these files. Also
upload both of these files to S3.

With that in place, prepare and upload the final configuration to S3 with::

    bcbio_vm.py template s3://your-project/your-analysis/template.yaml s3://your-project/your-analysis/name.csv

This will find the input files in the ``s3://your-project/your-analysis`` bucket, associate
fastq and BAM files with the right samples, and add a found BED files as
``variant_regions`` in the configuration. It will then upload the final
configuration back to S3 as ``s3://your-project/your-analysis/name.yaml``, which you can run
directly from a bcbio cluster on AWS. By default, bcbio will use the us-east S3
region, but you can specify a different region in the s3 path to the
metadata file: ``s3://your-project@eu-central-1/your-analysis/name.csv``

We currently support human analysis with both the GRCh37 and hg19 genomes. We
can also add additional genomes as needed by the community and generally welcome
feedback and comments on reference data support.

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

  bcbio_vm.py aws iam --region=us-east-1

The second configures a VPC to host bcbio::

  bcbio_vm.py aws vpc --region=us-east-1

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

The cluster will take five to ten minutes to start and be provisioned. If you encounter any
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

Create your project directory and link the global bcbio configuration file in there with:

- NFS file system (no Lustre)::

    mkdir /encrypted/your-project
    cd !$ && mkdir work && cd work

- Lustre file system::

    sudo mkdir /scratch/cancer-dream-syn3-exome
    sudo chown ubuntu !$
    cd !$ && mkdir work && cd work

If you started a single machine, run with::

    bcbio_vm.py run -n 8 s3://your-project/your-analysis/name.yaml

Where the ``-n`` argument should be the number of cores on the machine.

To run on a full cluster::

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

By default the collectl stats will be in ``monitoring/collectl`` and plots in
``monitoring/graphs`` based on the above log timeframe. If you need to re-run
plots later after shutting the cluster down, you can use the `none` cluster flag
by running ``bcbio_vm.py graph bcbio-nextgen.log --cluster none``.

If you'd like to run graphing from a local non-AWS run, such as a local HPC cluster,
run ``bcbio_vm.py graph bcbio-nextgen.log --cluster local`` instead.

For convenience, there's a "serialize" flag ('-s') that saves the dataframe used
for plotting. In order to explore the data and extract specific datapoints
or zoom, one could just deserialize the output like a python pickle file:

```
    import cPickle as pickle
    with gzip.open("./monitoring/collectl_info.pickle.gz", "rb") as decomp:
        collectl_info = pickle.load(decomp)
        data, hardware, steps = collectl_info[1][0], collectl_info[1][1], collectl_info[1][2]
```

And plot, slice, zoom it in an jupyter notebook using matplotlib,
[highcharts](https://github.com/arnoutaertgeerts/python-highcharts).

In addition to plots, the
`summarize_timing.py <https://github.com/bcbio/bcbio-nextgen/blob/master/scripts/utils/summarize_timing.py>`_
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
also see the `latest example configuration <https://github.com/bcbio/bcbio-nextgen-vm/blob/master/elasticluster/config>`_.
in the bcbio-vm GitHub repository for more details on the other available options.
