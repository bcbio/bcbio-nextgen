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
to reduce costs and use spot or preemptible instances at the cost of requiring
more hands on attention during runs. If runs fail or need custom downstream
analysis you can stop the larger run instances and work with a smaller, cheaper
machine.

## Cloud provider setup

### Amazon Web Services

Tools used on your local machine:

- [Ansible](http://docs.ansible.com/ansible/intro_installation.html) with
  [Dependencies and environmental variables for AWS access](http://docs.ansible.com/ansible/guide_aws.html)
  -- automate starting up instances
- [saws](https://github.com/donnemartin/saws) -- manage and query running
  instances.

You can install these into an isolated conda environment and setup with:

    wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
    bash Miniconda2-latest-Linux-x86_64.sh -b -p tools
    ./tools/bin/pip install ansible saws
    ./tools/bin/aws configure

You'll need to create some basic AWS infrastructure to do runs. You can use the
automation in
[bcbio-vm](http://bcbio-nextgen.readthedocs.io/en/latest/contents/cloud.html#aws-setup)
to create these:

    bcbio_vm.py aws iam
    bcbio_vm.py aws vpc

or use the [AWS console](https://aws.amazon.com/) or saws. You need:

- An AWS Virtual Private Cloud (VPC). A default VPC is fine.
- A security group allowing port 22 ssh access to the machines.
- The name of a keypair to use for ssh access, where you have the private key
  stored locally.
- Optionally, an IAM role that allows access to S3 resources. This makes it
  easier to push/pull data to the instance.

Finally, create a volume that will contain the run and bcbio installation. It
should be in the same availability zone as your created VPC. You can create this
in the AWS console in the Volumes tab, or using saws/aws:

       aws ec2 create-volume --size 300 --availability-zone us-east-1d --encrypted
       aws ec2 create-tags --resources vol-00df42a6 --tags Key=Name,Value=exome-validation

Use this information to create a configuration file called `project_vars.yaml`:

    instance_type: t2.small
    spot_price: null
    image_id: ami-2d39803a
    vpc_subnet: subnet-0817ce51
    volume: vol-79ce2ede
    security_group: bcbio_cluster_sg
    keypair: kunkel-keypair
    iam_role: bcbio_full_s3_access
    region: us-east-1

With this in place you can launch your instance with:

    ansible-playbook -vvv launch_aws.yaml

This creates the instance, attaches the data volume, mounts the volume as
`/mnt/work` and installs basic system tools for running. On the first run for a
new data volume, this will not work cleanly since the filesystem is not prepared
and can't be mounted. The machine will be setup and you should ssh in and follow
the instructions below in 'Running bcbio' to create an ext4 filesystem on the
attached disk. Subsequent restarts with the same attached disk will then work
without any manual steps.

Get the Public DNS name of the created machine with:

    saws> aws ec2 describe-instances | grep Public

Then ssh into your machine with:

    ssh -i /path/to/your/private.key ubuntu@ec2-XX-XXX-XXX-XXX.compute-1.amazonaws.com

You can skip ahead to the section about running bcbio, which is the same across
all cloud providers. When finished, terminate the current instance with:

    saws> aws ec2 terminate-instances --instance-ids i-xxxxxxx

Your bcbio installation and analysis is in the separate data volume. When
finished you can snapshot this for long term storage.

### Google Compute

- [Ansible](http://docs.ansible.com/ansible/intro_installation.html) with
  [dependencies and environmental variables for Google Compute access](http://docs.ansible.com/ansible/guide_gce.html)
  -- automate starting up instances
- [gloud from the Google Cloud SDK](https://cloud.google.com/sdk/) -- command
  line interface to access and manage instances. You can install with
  `bcbio_conda install -c bioconda google-cloud-sdk`

[Console](https://console.cloud.google.com)

Launch your instance with:

    ansible-playbook -vvv launch_gce.yaml

Then access the machine using your run name:

    gcloud compute ssh ubuntu@giab-val-work

When finished, you can terminate the instance with:

    gcloud compute instances delete giab-val-work

### Microsoft Azure

Tools used on your local machine:

- [Ansible](http://docs.ansible.com/ansible/intro_installation.html) with
  [dependencies and service principal based environmental variables for Azure access](https://docs.ansible.com/ansible/guide_azure.html).
  -- automate starting up instances
- [azure cross platform command line interface](https://github.com/Azure/azure-xplat-cli#features)
  -- command line interface to access and manage instances. You can install with
  `bcbio_conda install -c bioconda azure-cli`

Create a file with your instance and run information called `project_vars.yaml`:

    instance_type: Standard_DS1_v2
    run_name: giab-val-work
    storage_account_name: giabvalidation
    resource_group: giab-validation
    volume: giab-val-data
    image_id:
      publisher: Canonical
      offer: UbuntuServer
      sku: 16.04.0-LTS
      version: latest
    ssh_public_keys:
      - path: '/home/ubuntu/.ssh/authorized_keys'
        key_data: 'ssh-rsa Add your public key for connecting here'

Launch your instance with:

    ansible-playbook -vvv launch_azure.yaml

Then access the machine by finding the IP address and using ssh to attach:

    azure vm show giab-validation giab-val-work | grep Public
    ssh -i /path/to/your/private.key ubuntu@52.186.126.145

When finished, you can terminate the instance with:

    azure vm delete giab-validation giab-val-work

## Running bcbio

On the first run you'll need to create a project directory to work in:

    sudo mkdir /mnt/work/your-project
    sudo chown ubuntu /mnt/work/your-project

and [install bcbio](http://bcbio-nextgen.readthedocs.io/en/latest/contents/installation.html)
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
