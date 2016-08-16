# Simplified bcbio cloud usage

Many bcbio analyses can run on a single machine. A powerful multicore
instance can run smaller numbers of samples, targeted sequencing, or non-time
critical processing tasks without the overhead of spinning up a cluster and
shared filesystem.

This guide describes a simplified approach to run bcbio on
[Amazon Web Services](https://aws.amazon.com/),
[Google Cloud](https://cloud.google.com/) and
[Microsoft Azure](https://azure.microsoft.com):

- Install bcbio and your analysis on a persistent data volume.
- Use an ansible script to launch an instance with the data volume attached.
- Run bcbio analyses on the instance, stopping and changing instance sizes as
  needed to complete the run.

This assumes familiarity with the cloud platform you choose and offers minimal
automation compared to other
[Cloud integration](http://bcbio-nextgen.readthedocs.io/en/latest/contents/cloud.html)
within bcbio. The additional control over managing resources allows researchers
to reduce costs and use spot or interruptable instances at the cost of requiring
more hands on attention during runs. If runs fail or need custom downstream
analysis you can stop the larger run instances and work with a smaller, cheaper
machine.

## Cloud provider setup

### Amazon Web Services

Tools used on your local machine:

- [Ansible](http://docs.ansible.com/ansible/intro_installation.html) with
  [Dependencies and environmental variables set](http://docs.ansible.com/ansible/guide_aws.html)
  -- automate starting up instances
- [saws](https://github.com/donnemartin/saws) -- manage and query running instances

You'll need to use the [AWS console](https://aws.amazon.com/) or saws to setup
some basic infrastructure to do runs:

- An AWS Virtual Private Cloud (VPC). A default VPC is fine.
- A security group allowing port 22 ssh access to the machines.
- The name of a keypair to use for ssh access, where you have the private key
  stored locally.
- A volume that will contain the run and bcbio installation.
- Optionally, an IAM role that allows access to S3 resources. This makes it
  easier to push/pull data to the instance.

Use this information to create a configuration file called `project_vars.yaml`:

    instance_type: t2.small
    spot_price: null
    image_id: ami-2d39803a
    vpc_subnet: subnet-0817ce51
    volume: vol-79ce2ede
    security_group: bcbio_cluster_sg
    keypair: kunkel-keypair
    iam_role: bcbio_full_s3_access

With this in place you can launch your instance with:

    ansible-playbook -vvv launch_aws.yaml

This creates the instance, attaches the data volume, mounts the volume as
`/mnt/work` and installs basic system tools for running. Get the Public DNS
name of the created machine with:

    saws> aws ec2 describe-instances | grep Public

Then ssh into your machine with:

    ssh -i /path/to/your/private.key ubuntu@ec2-XX-XXX-XXX-XXX.compute-1.amazonaws.com

You can skip ahead to the section about running bcbio, which is the same across
all cloud providers. When finished, terminate the current instance with:

    saws> aws ec2 terminate-instances --instance-ids i-xxxxxxx

Your bcbio installation and analysis is in the separate data volume. When
finished you can snapshot this for long term storage.

### Google Compute

### Microsoft Azure

## Running bcbio

On the first run you'll need to create a filesystem on the data volume and setup
a directory for bcbio and your project:

    sudo mkfs -t ext4 /dev/xvdf
    sudo mkdir /mnt/work/bcbio
    sudo chown ubuntu /mnt/work/bcbio
    sudo mkdir /mnt/work/your-project
    sudo chown ubuntu /mnt/work/your-project

Then [install bcbio](http://bcbio-nextgen.readthedocs.io/en/latest/contents/installation.html)
on the working volume with the genomes and aligner indices you need:

    wget https://raw.github.com/chapmanb/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
    python bcbio_nextgen_install.py /mnt/work/bcbio --tooldir=/mnt/work/bcbio --genomes GRCh37 --aligners bwa

And you're ready to do an analysis in `/mnt/work/your-project`. Add your
samples, create a
[project configuration](http://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html#automated-sample-configuration)
and then run the analysis from a work directory.

Practically you can do all the bcbio installation and project setup with a
smaller instance, shutdown that instance and spin up a larger one that matches
your project needs for the actual run.

For runs, the ansible script adds the bcbio installation to your path so start
an analysis with:

    cd /mnt/work/your-project/work
    bcbio_nextgen.py ../config/your-project.yaml -n 16

To optimize resource usage, especially for single sample projects, edit
`/mnt/work/bcbio/galaxy/bcbio_system.yaml` so cores and memory match your instance:

    default:
      memory: 3G
      cores: 16
      jvm_opts: ["-Xms750m", "-Xmx3500m"]
