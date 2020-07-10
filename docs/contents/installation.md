# Installation

## Fresh installation (HPC cluster, server, AMI instance)

### 1. Install bcbio package and tools

`bcbio_nextgen_install.py` script installs:
- bcbio-nextgen python package;
- python library dependencies;
- third party analysis tools:

```shell
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
python bcbio_nextgen_install.py [bcbio_path] --tooldir=[bcbio_tools_path] --nodata 
```

You have to specify where to install bcbio in your filesystem and where to install tools, for example:
```
python bcbio_nextgen_install.py /bcbio --tooldir=/bcbio/tools --nodata 
```

or inside your home directory:
```
python bcbio_nextgen_install.py /home/user/bcbio --tooldir=/home/user/bcbio/tools --nodata 
```

Installation takes 30 minutes or more (depending on the speed of your storage and Internet connection). Recommended HPC job parameters for the installation process: 1 CPU core, 2GB memory, and 1 hour run time.

Check if installation works:
```
which bcbio_nextgen.py
bcbio_nextgen.py --version
```

### 2. Install data 

Bcbio needs reference files, indices, and databases to run any analyses. It is possible to install bcbio package and data at once, but we recommend to split these steps, because (i) some datatargets (dbNSFP, gnomad, snpEff) may take tens of hours or several days to finish, they could break in the middle due to unstable connections, i.e. it is better to tackle them one by one; (ii) you can re-use your data installation between bcbio instances. Data does not change much even between years, so you can just link `ln -s /old_bcbio/genomes /new_bcbio/genomes`

```
bcbio_nextgen.py -u skip --genomes hg38 --aligners bwa
```

This command installs hg38 human reference genome and bwa aligner index - the bare minimum required to run germline or somatic variant calling pipelines.

## Installation notes
- bcbio should install cleanly on Linux systems. For Mac OSX, we suggest trying [bcbio-vm](https://github.com/bcbio/bcbio-nextgen-vm) which runs bcbio on [Cloud](cloud) or isolates all the third party tools inside a Docker container. bcbio-vm is still a work in progress but not all of the dependencies bcbio uses install cleanly on OSX.
- Don't run the installer with sudo or as the root user. Do not use directories with `:` in the name, it is not POSIX compliant and will cause installation failures.
- To use custom mirrors for `conda-forge` and `bioconda` channels used during bcbio installation, set appropriate [channel alias](https://docs.conda.io/projects/conda/en/latest/user-guide/configuration/use-condarc.html#set-a-channel-alias-channel-alias) in your `.condarc` configuration file.
- The machine will need to have some basic requirements for installing and running bcbio:
  * Python 2.7 or Python 3.x
  * Basic system setup for unpacking files: tar, gzip, unzip, bzip2, xz-utils
  * The git version control system (<https://git-scm.com/>)
  * wget for file retrieval (<https://www.gnu.org/software/wget/>)
- Optional tool specific requirements:
  * Java 1.7, needed when running GATK < 3.6 or MuTect. This must be available in your path so typing `java -version`   resolves a 1.7 version. bcbio distributes Java 8 as part of the Anaconda installation for recent versions of GATK and MuTect2. You can override the Java 8 installed with bcbio by setting `BCBIO_JAVA_HOME=/path/to/your/javadir` if you have the Java you want in `/path/to/your/javadir/bin/java`.
  * An OpenGL library, like [Mesa](https://www.mesa3d.org/) (On Ubuntu/deb systems: `libglu1-mesa`, On RedHat/rpm systems: `mesa-libGLU-devel`). This is only required for cancer heterogeneity analysis with BubbleTree.
  * The Pisces tumor-only variant callers requires the [Microsoft .NET runtime](https://docs.microsoft.com/en-us/dotnet/core/install/linux-package-managers).
- The [bcbio-nextgen Dockerfile](https://github.com/bcbio/bcbio-nextgen/blob/master/Dockerfile#L5) contains the packages needed to install on bare Ubuntu systems.
- The automated installer creates a fully integrated environment that allows simultaneous updates of the framework, third party tools and biological data. This offers the advantage over manual installation of being able to manage and evolve a consistent analysis environment as algorithms continue to evolve and improve. Installing this way is as isolated and self-contained as possible without virtual machines or lightweight system containers like [Docker](https://www.docker.com/). 

## Installation parameters
Run 
```
bcbio_nextgen.py upgrade --help
```
to see all supported installation options:
```
bcbio_nextgen.py upgrade --help
usage: bcbio_nextgen.py upgrade [-h] [--cores CORES] [--tooldir TOOLDIR]
                                [--tools]
                                [-u {stable,development,system,deps,skip}]
                                [--toolconf TOOLCONF] [--revision REVISION]
                                [--toolplus TOOLPLUS]
                                [--datatarget {variation,rnaseq,smallrna,gemini,vep,dbnsfp,dbscsnv,battenberg,kraken,ericscript,gnomad}]
                                [--genomes {GRCh37,hg19,hg38,hg38-noalt,mm10,mm9,rn6,rn5,canFam3,dm3,galGal4,phix,pseudomonas_aeruginosa_ucbpp_pa14,sacCer3,TAIR10,WBcel235,xenTro3,GRCz10,GRCz11,Sscrofa11.1,BDGP6}]
                                [--aligners {bwa,rtg,hisat2,bbmap,bowtie,bowtie2,minimap2,novoalign,twobit,bismark,snap,star,seq}]
                                [--data] [--cwl] [--isolate]
                                [--distribution {ubuntu,debian,centos,scientificlinux,macosx}]

optional arguments:
  -h, --help            show this help message and exit
  --cores CORES         Number of cores to use if local indexing is necessary.
  --tooldir TOOLDIR     Directory to install 3rd party software tools. Leave
                        unspecified for no tools
  --tools               Boolean argument specifying upgrade of tools. Uses
                        previously saved install directory
  -u {stable,development,system,deps,skip}, --upgrade {stable,development,system,deps,skip}
                        Code version to upgrade
  --toolconf TOOLCONF   YAML configuration file of tools to install
  --revision REVISION   Specify a git commit hash or tag to install
  --toolplus TOOLPLUS   Specify additional tool categories to install
  --datatarget {variation,rnaseq,smallrna,gemini,vep,dbnsfp,dbscsnv,battenberg,kraken,ericscript,gnomad}
                        Data to install. Allows customization or install of
                        extra data.
  --genomes {GRCh37,hg19,hg38,hg38-noalt,mm10,mm9,rn6,rn5,canFam3,dm3,galGal4,phix,pseudomonas_aeruginosa_ucbpp_pa14,sacCer3,TAIR10,WBcel235,xenTro3,GRCz10,GRCz11,Sscrofa11.1,BDGP6}
                        Genomes to download
  --aligners {bwa,rtg,hisat2,bbmap,bowtie,bowtie2,minimap2,novoalign,twobit,bismark,snap,star,seq}
                        Aligner indexes to download
  --data                Upgrade data dependencies
  --cwl                 Install code and data for running CWL workflows
  --isolate             Created an isolated installation without PATH updates
  --distribution {ubuntu,debian,centos,scientificlinux,macosx}
                        Operating system distribution
```

Some useful arguments are:
* `--isolate` Avoid updating the user's `~/.bashrc` if installing in a non-standard PATH. This facilitates creation of isolated modules without disrupting the user's environmental setup. Manually edit your `~/.bashrc` to allow bcbio runs with:
    ```shell
    export PATH=/path_to_bcbio/anaconda/bin:/path_to_bcbio/tools/bin:$PATH
    ```
* `--nodata` Do not install genome data.

## On a Virtual Machine

If you are looking to quickly try out bcbio-nextgen on your personal machine before installing it on your cluster, installing bcbio-nextgen on a virtual machine is easy using [Vagrant](https://www.vagrantup.com/).

### macOS

* Install [Git](https://git-scm.com/download/mac), [VirtualBox](https://download.virtualbox.org/virtualbox/6.1.6/VirtualBox-6.1.6-137129-OSX.dmg), and [Vagrant](https://releases.hashicorp.com/vagrant/2.2.7/vagrant_2.2.7_x86_64.dmg)
* Download bcbio-nextgen and provision Vagrant VM:
    ```shell
    git clone git@github.com:bcbio/bcbio-nextgen.git
    cd bcbio-nextgen
    vagrant up
    ```
* Install bcbio-nextgen (this should take about 30 minutes):
    ```shell
    vagrant ssh
    python /vagrant/scripts/bcbio_nextgen_install.py ~/local/share/bcbio --tooldir=~/local --nodata
    ```
Optional steps:
* Inside the VM (`vagrant ssh`):
  * Test your installation once it's complete:
    ```shell
    bcbio_nextgen.py --version
    ```
  * Set the time zone in the VM for easier log viewing, for example:
    ```shell
    sudo timedatectl set-timezone America/New_York
    ```
* Outside the VM:
  * To make any additional data from the host available inside the VM (for example: reference genomes, pipeline inputs, etc) set `BCBIO_DATA_DIR` environment variable on the host to a directory that contains the data, for example:
    ```shell
    export BCBIO_DATA_DIR=~/biodata
    vagrant reload
    ```
    This directory will be mounted inside Vagrant VM under `/data` 

## Upgrade

We use the same automated installation process for performing upgrades of tools, software and data in place. Since there are multiple targets and we want to avoid upgrading anything unexpectedly, we have specific arguments for each. Generally, you'd want to upgrade the code, tools and data together with:
```shell
bcbio_nextgen.py upgrade -u stable --tools --data
```

**2020-05-21: in bcbio 1.2.3 upgrade -u stable is broken, use -u development, or -u skip, this will be fixed in bcbio 1.2.4**

Tune the upgrade with these options:
* `-u` Type of upgrade to do for bcbio-nextgen code. `stable` gets the most recent released version and `development` retrieves the latest code from GitHub.
* `--datatarget` Customized installed data or download additional files not included by default: 
* `--toolplus` Specify additional tools to include. See the section on `extra software` for more details.
* `--genomes` and `--aligners` options add additional aligner indexes to download and prepare. `bcbio_nextgen.py upgrade -h` lists available genomes and aligners. If you want to install multiple genomes or aligners at once, specify `--genomes` or `--aligners` multiple times, like this: `--genomes GRCh37 --genomes mm10 --aligners bwa --aligners bowtie2`
* Leave out the `--tools` option if you don't want to upgrade third party tools. If using `--tools`, it will use the same directory as specified during installation. If you're using an older version that has not yet gone through a successful upgrade or installation and saved the tool directory, you should manually specify `--tooldir` for the first upgrade. You can also pass `--tooldir` to install to a different directory.
* Leave out the `--data` option if you don't want to get any upgrades of associated genome data.
* Some aligners such as STAR don't have pre-built indices due to the large file sizes of these. You set the number of cores to use for indexing with `--cores 8`.
* For example, recommended HPC job parameters for `bcbio_nextgen.py upgrade -u skip --data --datatarget rnaseq --genomes GRCh37` are: 2 CPU cores, 2GB memory, and 2 hours run time.

## Customizing data installation

bcbio supports the following [genome references](https://github.com/chapmanb/cloudbiolinux/blob/master/config/biodata.yaml), 12 of them have [additional data downloads](https://github.com/chapmanb/cloudbiolinux/tree/master/ggd-recipes). If you need a reference which is absent in the list, you may install it as a [custom genome](configuration.html#adding-custom-genomes). 

bcbio installs associated data files for sequence processing, and you're able to customize this to install larger files or change the defaults. Use the `--datatarget` flag (potentially multiple times) to customize or add new targets.

By default, bcbio will install data files for `variation`, `rnaseq` and `smallrna` but you can sub-select a single one of these if you don't require other analyses. The available targets are:
* `variation` -- Data files required for variant calling: SNPs, indels and structural variants. These include files for annotation like dbSNP, associated files for variant filtering, coverage and annotation files.
* `rnaseq` -- Transcripts and indices for running RNA-seq. The transcript files are also used for annotating and prioritizing structural variants.
* `smallrna` -- Data files for doing small RNA analysis.
* `gemini` -- The [GEMINI](https://gemini.readthedocs.io) framework associates publicly available metadata with called variants, and provides utilities for query and analysis. This target installs the required GEMINI data files, including [ExAC](http://exac.broadinstitute.org/).
* `gnomad` -- [gnomAD](https://gnomad.broadinstitute.org/) is a large scale collection of genome variants, expanding on ExAC to include whole genome and more exome inputs. This is a large 25Gb download, available for human genome builds GRCh37, hg19 and hg38.
* `vep` -- Data files for the [Variant Effects Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html). To use VEP as an alternative to the default installed snpEff, set `vep` in the `variant calling` configuration.
* `dbnsfp` -- Like CADD, [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP) provides integrated and generalized metrics from multiple sources to help with prioritizing variations for follow up. The files are large: dbNSFP is 10Gb, expanding to 100Gb during preparation.
* `dbscsnv` -- [dbscSNV](https://sites.google.com/site/jpopgen/dbNSFP) includes all potential human SNVs within splicing consensus regions (−3 to +8 at the 5' splice site and −12 to +2 at the 3' splice site), i.e. scSNVs, related functional annotations and two ensemble prediction scores for predicting their potential of altering splicing.
* `battenberg` -- Data files for [Battenberg](https://github.com/cancerit/cgpBattenberg), which detects subclonality and copy number changes in whole genome cancer samples.
* `kraken` -- Database for [Kraken](https://ccb.jhu.edu/software/kraken/), optionally used for contamination detection.
* `ericscript` -- Database for [EricScript](https://sites.google.com/site/bioericscript/), based gene fusion detection. Supports hg38, hg19 and GRCh37.

For somatic analyses, bcbio includes [COSMIC](https://cancer.sanger.ac.uk/cosmic) v68 for hg19 and GRCh37 only. Due to license restrictions, we cannot include updated versions of
this dataset and hg38 support with the installer. To prepare these datasets yourself you can use [a utility script shipped with cloudbiolinux](https://github.com/chapmanb/cloudbiolinux/blob/master/utils/prepare_cosmic.py) that downloads, sorts and merges the VCFs, then copies into your bcbio installation:
```shell
export COSMIC_USER="you@example.org"
export COSMIC_PASS="your_cosmic_password"
bcbio_python prepare_cosmic.py 89 /path/to/bcbio
```
`/path/to/bcbio/` here is the directory one up from the `genomes` directory.

## Extra software

We're not able to automatically install some useful tools due to licensing restrictions, so we provide a mechanism to manually download and add these to bcbio-nextgen during an upgrade with the `--toolplus` command line option.

### GATK and MuTect/MuTect2

bcbio includes an installation of GATK4, which is freely available for all uses. This is the default runner for HaplotypeCaller or MuTect2. If you want to use an older version of GATK, it requires manual installation. This is freely available for academic users, but requires a [license for commercial use](https://gatk.broadinstitute.org/hc/en-us#licensing). It is not freely redistributable, so requires a manual download from the [GATK download](https://console.cloud.google.com/storage/browser/gatk-software/package-archive) site, direct [link](https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2). You also need to include `tools_off: [gatk4]` in your configuration for runs: see `changing bcbio defaults`.

To install GATK3, register with the pre-installed gatk bioconda wrapper:
```shell
gatk3-register /path/to/GenomeAnalysisTK.tar.bz2
```
If you're not using the most recent post-3.6 version of GATK, or using a nightly build, you can add `--noversioncheck` to the command line to skip comparisons to the GATK version.

[MuTect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) is distributed with GATK in versions 3.5 and later.

To install versions of GATK < 3.6, download and unzip the latest version from the GATK distribution. Then make this jar available to bcbio-nextgen with:
```shell
bcbio_nextgen.py upgrade --tools --toolplus gatk=/path/to/gatk/GenomeAnalysisTK.jar
```
This will copy the jar and update your `bcbio_system.yaml` and manifest files to reflect the new version.

MuTect also has similar licensing terms and requires a license for commercial use. After [downloading the MuTect jar](https://www.broadinstitute.org/gatk/download/), make it available to bcbio:
```shell
bcbio_nextgen.py upgrade --tools --toolplus mutect=/path/to/mutect/mutect-1.1.7.jar
```
Note that muTect does not provide an easy way to query for the current version, so your input jar needs to include the version in the name.

## System requirements

bcbio-nextgen provides a wrapper around external tools and data, so the actual tools used drive the system requirements. For small projects, it should install on workstations or laptops with a couple GB of memory, and then scale as needed on clusters or multicore machines.

Disk space requirement for the tools, including all system packages is about 22GB (or more, depending on the type of the file system). Biological data requirements will depend on the genomes and aligner indices used, but a suggested install with GRCh37 and bowtie/bwa2 indexes uses approximately 35GB of storage during preparation and ~25GB after:
```shell
$ du -shc genomes/Hsapiens/GRCh37/*
3.8G  bowtie2
5.1G  bwa
3.0G  rnaseq-2014-05-02
3.0G  seq
340M  snpeff
4.2G  variation
4.4G  vep
23.5G total
```

## Troubleshooting

### Proxy or firewall problems

Some steps retrieve third party tools from GitHub, which can run into issues if you're behind a proxy or block git ports. To instruct git to use `https://` globally instead of `git://`:
```shell
git config --global url.https://github.com/.insteadOf git://github.com/
```

### GATK or Java Errors

Most software tools used by bcbio require Java 1.8. bcbio distributes an OpenJDK Java build and uses it so you don't need to install anything. Older versions of GATK (< 3.6) and MuTect require a locally installed Java 1.7. If you have version incompatibilities, you'll see errors like:
```
Unsupported major.minor version 51.0
```
Fixing this requires either installing Java 1.7 for old GATK and MuTect or avoiding pointing to an incorrect java (`unset JAVA_HOME`). You can also tweak the java used by bcbio, described in the [Automated](#automated) installation section.

### ImportErrors

Import errors with tracebacks containing Python libraries outside of the bcbio distribution (`/path/to/bcbio/anaconda`) are often due to other conflicting Python installations. bcbio tries to isolate itself as much as possible but external libraries can get included during installation due to the `PYTHONHOME` or `PYTHONPATH` environmental variables or local site libraries. These commands will temporary unset those to get bcbio installed, after which it should ignore them automatically:
```shell
unset PYTHONHOME
unset PYTHONPATH
export PYTHONNOUSERSITE=1
```
Finally, having a `.pydistutils.cfg` file in your home directory can mess with where the libraries get installed. If you have this file in your home directory, temporarily renaming it to something else may fix your installation issue.

## Manual process

The manual process does not allow the in-place updates and management of third party tools that the automated installer makes possible. It's a more error-prone and labor intensive process. If you find you can't use the installer we'd love to hear why to make it more amenable to your system. If you'd like to develop against a bcbio installation, see the documentation on setting up a `development environment`.

### Tool requirements

The code drives a number of next-generation sequencing analysis tools that you need to install on any machines involved in the processing. The [CloudBioLinux](http://cloudbiolinux.org) toolkit provides automated scripts to help with installation for both software and associated data files:
```shell
fab -f cloudbiolinux/fabfile.py -H localhost install_biolinux:flavor=ngs_pipeline_minimal
```
You can also install them manually, adjusting locations in the `resources` section of your `bcbio_system.yaml` configuration file as needed. The CloudBioLinux infrastructure provides a full list of third party software installed with bcbio-nextgen in [packages-conda.yaml](https://github.com/chapmanb/cloudbiolinux/blob/master/contrib/flavor/ngs_pipeline_minimal/packages-conda.yaml), which lists all third party tools installed through [Bioconda](https://bioconda.github.io/).

### Data requirements

In addition to existing bioinformatics software the pipeline requires associated data files for reference genomes, including pre-built indexes for aligners. The [CloudBioLinux](http://cloudbiolinux.org) toolkit again provides an automated way to download and prepare these reference genomes:
```shell
fab -f data_fabfile.py -H localhost -c your_fabricrc.txt install_data_s3:your_biodata.yaml
```
The [biodata.yaml](https://github.com/chapmanb/cloudbiolinux/blob/master/config/biodata.yaml) file contains information about what genomes to download. The [fabricrc.txt](https://github.com/chapmanb/cloudbiolinux/blob/master/config/fabricrc.txt) describes where to install the genomes by adjusting the `data_files` variable. This creates a tree structure that includes a set of Galaxy-style location files to describe locations of indexes:
```
├── galaxy
│   ├── tool-data
│   │   ├── alignseq.loc
│   │   ├── bowtie_indices.loc
│   │   ├── bwa_index.loc
│   │   ├── sam_fa_indices.loc
│   │   └── twobit.loc
│   └── tool_data_table_conf.xml
├── genomes
│   ├── Hsapiens
│   │   ├── GRCh37
│   │   └── hg19
│   └── phiX174
│       └── phix
└── liftOver
```
Individual genome directories contain indexes for aligners in individual sub-directories prefixed by the aligner name. This structured scheme helps manage aligners that don't have native Galaxy `.loc` files. The automated installer will download and set this up
automatically:
```
`-- phix
    |-- bowtie
    |   |-- phix.1.ebwt
    |   |-- phix.2.ebwt
    |   |-- phix.3.ebwt
    |   |-- phix.4.ebwt
    |   |-- phix.rev.1.ebwt
    |   `-- phix.rev.2.ebwt
    |-- bowtie2
    |   |-- phix.1.bt2
    |   |-- phix.2.bt2
    |   |-- phix.3.bt2
    |   |-- phix.4.bt2
    |   |-- phix.rev.1.bt2
    |   `-- phix.rev.2.bt2
    |-- bwa
    |   |-- phix.fa.amb
    |   |-- phix.fa.ann
    |   |-- phix.fa.bwt
    |   |-- phix.fa.pac
    |   |-- phix.fa.rbwt
    |   |-- phix.fa.rpac
    |   |-- phix.fa.rsa
    |   `-- phix.fa.sa
    |-- novoalign
    |   `-- phix
    |-- seq
    |   |-- phix.dict
    |   |-- phix.fa
    |   `-- phix.fa.fai
    `-- ucsc
        `-- phix.2bit
```
