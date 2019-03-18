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
- [bcbio-vm](http://bcbio-nextgen.readthedocs.io/en/latest/contents/cloud.html#aws-setup)
  -- automates creation of AWS resources

Install these into an isolated conda environment and setup with:

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p tools
    ./tools/bin/conda install -c conda-forge -c bioconda python=3 bcbio-nextgen-vm
    ./tools/bin/pip install ansible saws boto
    ./tools/bin/aws configure

Provide AWS access for bcbio-vm, ansible and saws using [IAM to create (or find)
your access keys](https://docs.aws.amazon.com/IAM/latest/UserGuide/id_credentials_access-keys.html):


    export AWS_ACCESS_KEY_ID='AK123'
    export AWS_SECRET_ACCESS_KEY='abc123'

bcbio-vm has an automated script to setup the AWS infrastructure from running:

    bcbio_vm.py aws ansible us-east-1d --keypair

Replace `us-east-1d` with the AWS availability zone you'd like to run in. This
creates a `project_vars.yaml` file with the following information:

- An AWS image ID for your region, selected from
  [Ubuntu Amazon EC2 AMI Locator](http://cloud-images.ubuntu.com/locator/ec2/)
  using an option with Instance Type `hvm:ebs-ssd`.
- An AWS Virtual Private Cloud (VPC) subnet.
- A security group allowing port 22 ssh access to the machines.
- The name of a [keypair](https://console.aws.amazon.com/ec2/v2/home?region=us-east-1#KeyPairs:sort=keyName) 
  to use for ssh access, where you have the private key stored locally. If you
  used the bcbio-vm automated setup, you'll have a private keypair in
  `~/.bcbio/aws_keypairs/bcbio`.
- An IAM role that allows access to S3 resources. This makes it
  easier to push/pull data to the instance.

You can also use the [AWS console](https://aws.amazon.com/) or saws to create
these manually.

Create a volume that will contain the run and bcbio installation. It
should be in the same availability zone as specified above. You can create this
in the AWS console in the Volumes tab, or using saws/aws:

    aws ec2 create-volume --encrypted --volume-type gp2 --size 300 --availability-zone us-east-1d
    aws ec2 create-tags --resources vol-00df42a6 --tags Key=Name,Value=exome-validation

Add your volume ID to the `project_vars.yaml` configuration file:

    instance_type: t2.small
    spot_price: null
    volume: vol-50ed15c1
    keypair: bcbio
    image_id: ami-6edd3078
    vpc_subnet: subnet-3628576d
    iam_role: bcbio_full_s3_access
    security_group: bcbio_cluster_sg
    region: us-east-1
    zone: us-east-1d

Edit `~/.ssh/config` to enable clean connections to AWS instances:

    Host *.amazonaws.com
      StrictHostKeyChecking no
      UserKnownHostsFile=/dev/null
      IdentityFile /path/to/aws_keypairs/bcbio
      IdentitiesOnly yes
      ControlPath ~/.ssh/%r@%h:%p

With this in place you can launch your instance with:

    ansible-playbook -i 'localhost,' -vvv launch_aws.yaml

This creates the instance, attaches the data volume, mounts the volume as
`/mnt/work` and installs basic system tools. Get the Public DNS name of the
created machine with:

    saws> aws ec2 describe-instances | grep Public

Then ssh into your machine with:

    ssh ubuntu@ec2-XX-XXX-XXX-XXX.compute-1.amazonaws.com

On the first run with a new data volume, install bcbio and setup a work
directory for running following the instructions in 'Running bcbio.' Then you're
ready to run an analysis. This on-machine setup is the same across all cloud
providers. When finished, terminate the current instance with:

    saws> aws ec2 terminate-instances --instance-ids i-xxxxxxx

Your bcbio installation and analysis is in the separate data volume. When
finished you can snapshot this for long term storage.

### Google Compute

Tools used on your local machine:

- [Ansible](http://docs.ansible.com/ansible/intro_installation.html) with
  [dependencies and environmental variables for Google Compute access](http://docs.ansible.com/ansible/guide_gce.html)
  -- automate starting up instances
- [gloud from the Google Cloud SDK](https://cloud.google.com/sdk/) -- command
  line interface to access and manage instances. You can install with
  `bcbio_conda install -c conda-forge -c bioconda google-cloud-sdk`
- [Web based console](https://console.cloud.google.com)

Use the console to create a project to hold your analysis
, then add [ssh key access to your project](https://cloud.google.com/compute/docs/instances/adding-removing-ssh-keys#project-wide).
Then locally, [log into gcloud](https://cloud.google.com/compute/docs/gcloud-compute/)
and select your project:

    gcloud init

Create a disk to store bcbio and associated data:

    gcloud compute disks create dv-bcbio-vol --size 250GB --type pd-ssd --zone us-east1-b

Finally [create credentials](http://docs.ansible.com/ansible/latest/scenario_guides/guide_gce.html#credentials)
for connecting to your instance. These are key pairs used to automatically
authenticate to created instances. To create, go to "APIs and Services",
"Credentials", "Create credentials" and finally "Service account key", then
download the json file with key pairs after creating it.
You'll need the e-mail associated with the service account, found under "IAM and
Admin."

Use this information to create a `project_vars.yaml` configuration file:

    instance_type: n1-standard-1
    image_id: ubuntu-1604-xenial-v20180323
    zone: us-east1-b
    service_account_email: 129107966647-compute@developer.gserviceaccount.com
    credentials_file: /home/chapmanb/.ssh/gce/deepvariant-trustedtester-d4f7663f4adf.json
    project_id: deepvariant-trustedtester
    volume: dv-bcbio-vol
    run_name: dv-bcbio

Launch your instance with:

    ansible-playbook -vvv launch_gce.yaml

Then access the machine using your run name:

    gcloud compute ssh ubuntu@dv-bcbio

When finished, you can terminate the instance with:

    gcloud compute instances delete dv-bcbio

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

    wget https://raw.github.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
    python bcbio_nextgen_install.py /mnt/work/bcbio --tooldir=/mnt/work/bcbio --genomes GRCh37 --aligners bwa

To run CWL, you'll also want to [install
bcbio-vm](https://bcbio-nextgen.readthedocs.io/en/latest/contents/cwl.html#getting-started)
with:

    export TARGETDIR=/mnt/work/bcbio/bcbio-vm
    export BINDIR=/mnt/work/bcbio/bin

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

# Running with Toil

This in progress documentation describes running bcbio generated CWL using
[Toil autoscaling on AWS](http://toil.readthedocs.io/en/latest/running/cloud.html).

- Find a [Docker Toil tag](https://quay.io/repository/ucsc_cgl/toil?tag=latest&tab=tags)
  corresponding to the version you're running (if not a stable release):

        curl -s 'https://quay.io/api/v1/repository/ucsc_cgl/toil/tag/?limit=20' \
          | jq '.tags[] | (.name)'

- Pick a VPC subnet to run in. This subnet needs to enable auto-assign public IP
  addresses.

- Launch the cluster

        export TOIL_AWS_ZONE=us-east-1d
        export TOIL_APPLIANCE_SELF='quay.io/ucsc_cgl/toil:3.7.0a1.dev375-587871c5f4f219877af7d8e4f0a1c9544c510e65'
        toil launch-cluster -p aws bcbio --nodeType=t2.small --keyPairName=bcbio \
          --vpcSubnet subnet-3628576d

- Find local IP of machine (for Mesos LEADER_PRIVATE_IP):

        aws ec2 describe-instances | grep PrivateIp

- Attach to Docker container in head node:

        toil ssh-cluster toil-bcbio

- Run bcbio CWL toil test

        mkdir /home/run && cd /home/run
        wget -O test_bcbio_cwl.tar.gz https://github.com/bcbio/test_bcbio_cwl/archive/master.tar.gz
        tar -xzvpf test_bcbio_cwl.tar.gz
        cd test_bcbio_cwl-master
        <edit run_toil_aws.sh to change LEADER_PRIVATE_IP, from above) and JOB_STORE s3 bucket>
        bash run_toil_aws.sh

Current discussion topics:

  - Need to be able to specify root volume sizes. Current 50Gb hardcoded default
    will fail for larger WGS samples and some test samples.
  - How does file staging work for shared files used in multiple steps, like
    BAMs for variant calling? Are they staged once on an AWS machine and re-used
    or pulled down multiple times? I'm running into disk space errors during
    runs which indicate they get pulled down multiple times and potentially not
    removed, or we have two jobs on the same machine with larger files.
    Re-starting the pipeline allows these jobs to finish.
  - Restarting jobs with/without `--no-container` while reusing job store. In
    general, can we adjust Toil parameters without needing to re-setup the
    jobStore? 
  - Can we scale down proceses that request too many cores (8 cores on 4 core
    machine)?
  - Similarly, how does passing files back to the S3 file store work? The
    pipeline current stalls on `batch_for_variantcall` which isn't doing much
    processing work but does aggressively split to run variant calls in
    parallel, so it return a large number of outputs each of which has the BAM
    file. I don't see these written to the filestore on S3 but am trying to
    brainstorm reason why this step takes a very long time to run.
  - Consider how to deal with input files already present in external S3
    buckets. Could we avoid copying these into the Toil work bucket if in the
    right region?
  - cgcloud versioning, need development version of Toil: https://github.com/BD2KGenomics/toil/issues/1458
  - Swap to using Ubuntu instead of relying on manual CoreOS Amazon approvals.

## Debugging tips

- Look up an appropriate CoreOS AMI for
  [your zone](https://github.com/BD2KGenomics/toil/issues/1470). This is only
  necessary if Toil can't automatically find an image:

        aws ec2 describe-images \
          --output table --query 'Images[*].{AMI:ImageId,Description:Description}' \
          --filters Name=owner-id,Values=679593333241 Name=description,Values='*stable*'
        export TOIL_AWS_AMI=ami-61659e77
