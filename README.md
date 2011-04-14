## Overview

Python scripts and modules for automated next gen sequencing analysis.
These provide a fully automated pipeline for taking sequencing results
from an Illumina sequencer, converting them to standard Fastq format,
aligning to a reference genome, doing SNP calling, and producing a
summary PDF of results.

The scripts are tightly integrated with the [Galaxy][o1]
web-based analysis tool. Samples are entered and tracked through a LIMS
system and processed results are uploading into Galaxy Data Libraries for
researcher access and additional analysis. See the
[installation instructions for the front end][o2] and a
[detailed description of the full system][o3].

![Overview of workflow][o4]

[o1]: http://galaxy.psu.edu/
[o2]: https://bitbucket.org/galaxy/galaxy-central/wiki/LIMS/nglims
[o3]: http://bcbio.wordpress.com/2011/01/11/next-generation-sequencing-information-management-and-analysis-system-for-galaxy/
[o4]: http://chapmanb.github.com/bcbb/nglims_organization.png

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

### Required libraries and data files

The code drives a number of next-generation sequencing analysis tools,
which will need to be installed along with associated data files. The
[CloudBioLinux][i2] and [Galaxy CloudMan][i3] projects contain
automated scripts to help with installation:

* [Install bioinformatics software][i4]
* Install data files like genome builds in association with
  Galaxy: [my script][i5] and [from the Galaxy team][i6].

### Scripts and configuration

Clone a copy of the code from GitHub:

      git clone git://github.com/chapmanb/bcbb.git

Install using the standard python approach. You will need a recent
version of Python 2 (2.6 or 2.7). The required python library
dependencies will also be installed:

      cd bcbb/nextgen && python setup.py install

Copy the YAML & ini files in config and adjust them to match your
environment. It is also a good idea to set your $PATH pointing to
any third-party binaries you are using.

### Testing

The test suite exercises the various scripts driving the analysis, so
are a good starting point to ensure everything has been installed
correctly. Tests are run from the main directory using [nose][i7]:

      nosetest -v -s

`tests/test_automated_analysis.py` is quite extensive and exercises
the full framework using an automatically downloaded test dataset. It
runs through barcode deconvolution, alignment and full SNP
analysis. The configuration for the tests can be tweaked for your
environment:

* `tests/data/automated/post_process.yaml` -- May need adjustment to
  point to installed software in non-standard locations.
* `tests/data/automated/run_info.yaml` -- The `analysis` variable
  can be changed to 'Standard' if SNP calling is not required in your
  environment. This will run a smaller pipeline of alignment and analysis.

### RabbitMQ messaging server

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

### ssh keys

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

### Sequencing machines

The backend has been developed using Illumina GAII and HiSeq sequencing
machines. It is designed to be general and support other platforms, and 
we welcome feedback from researchers with different machines at their
institutions.

Illumina machines produce run directories that include the date, machine
identifier, and flowcell ID:

    110216_HWI-EAS264_00035_FC638NPAAXX

The data and flowcell ID will be parsed from the directory named used as a
shortened reference name that includes an easier to remember date and the unique
flowcell ID for traceability.

### Development environment

The installation instructions assume that you have full root access to install
python modules and packages (production environment). If this is not the case,
you may want to install a python VirtualEnv and other tools automatically on your $HOME
to ease your development needs using the following script:

        wget https://bitbucket.org/brainstorm/custom_env/raw/1cd4f4ae27d5/pyHost.sh && ./pyHost.sh


## Internals: files generated by this pipeline

### Initial Fastq files (pre-analysis)

After basecalling, a number of FastQ files are generated via <a href="https://github.com/brainstorm/bcbb/blob/master/nextgen/scripts/solexa_qseq_to_fastq.py">solexa_qseq_to_fastq.py</a> script:

<pre>
1_081227_B45GT6ABXX_1.fastq
1_081227_B45GT6ABXX_2.fastq
(...)
</pre>

The template of them being:

lane_date_fcid_<1|2>.fastq

Where 1|2 is the forward and reverse read respectively.

### Post-processing generated files

Those are distributed in three directories: alignments, _barcode and the top level run directory.

#### Top level directory

Those files comprise both plain text files, images and structured data that can be useful both automatically and in human-readable form to get different aspects about the sequencing results.

###### run_summary.yaml

Contains a structured view of the run global parameters.

<pre>
- barcode_id: '2'
  barcode_type: SampleSheet
  lane: '8'
  metrics:
    Aligned: '1204752'
    Pair duplicates: '881332'
    Read length: '101'
    Reads: '11594437'
  request: ''
  researcher: ''
  sample: ''
</pre>

###### PDF files

Taking a look at an specific sample in a lane, we find the different sub-components that conform the final summary report (*-summary.pdf). The summary contains GC content, read quality figures and an insert size histogram.

<pre>
6_110126_B816J0ABXX_5_1_fastq_txt_quality.pdf
6_110126_B816J0ABXX_5_2_fastq.txt.segments.hist.pdf
6_110126_B816J0ABXX_5-sort-dup-gc.pdf
6_110126_B816J0ABXX_5_1_fastq.txt.segments.hist.pdf
6_110126_B816J0ABXX_5-sort_1_fastq_qual.pdf
6_110126_B816J0ABXX_5-sort-dup-insert.pdf
6_110126_B816J0ABXX_5_2_fastq_txt_quality.pdf
6_110126_B816J0ABXX_5-sort_2_fastq_qual.pdf
6_110126_B816J0ABXX_5-sort-summary.pdf
</pre>

###### BigWig files

Per-sample wigtoBigWig converted files.

##### Derived files

###### *-sort_1_fastq_stats.txt files

<table>
<tr>
<th>column</th>	<th>count</th>	<th>min</th>	<th>max</th>	<th>sum</th>	<th>mean</th>	<th>Q1</th>	<th>med</th>	<th>Q3</th>	<th>IQR</th>	<th>lW</th>	<th>rW</th>	<th>A_Count</th>	<th>C_Count</th>	<th>G_Count</th>	<th>T_Count</th>	<th>N_Count</th>	<th>Max_count</th>
</tr>
<tr><td>1</td>	<td>17837107</td>	<td>2</td>	<td>40</td>	<td>659805683</td>	<td>36.99</td>	<td>37</td>	<td>39</td>	<td>39</td>	<td>2</td>	<td>34</td>	<td>40</td>	<td>5303108</td>	<td>5298743</td>	<td>6299525</td>	<td>932665</td>	<td>3066</td>	<td>17837107</td></tr>
<tr><td>2</td>	<td>17837107</td>	<td>2</td>	<td>40</td>	<td>644870452	36</td><td>.15</td><td>38</td>	<td>39</td>	<td>3</td>	<td>32</td>	<td>40</td>	<td>12461524</td>	<td>1264534</td>	<td>3067874</td>	<td>1043174</td>	<td>1</td>	<td>17837107</td></tr>
</table>

###### *_metrics files

They contain plain text values from picard, namely: 

<pre>
1_110126_B816J0ABXX_5-sort-dup.align_metrics
1_110126_B816J0ABXX_5-sort-dup.gc_metrics
1_110126_B816J0ABXX_5-sort-dup.dup_metrics
1_110126_B816J0ABXX_5-sort-dup.insert_metrics
</pre>

They contain output metrics from Picard, which are parsed later and used to generate plots. For instance, for "sort-dup.align_metrics", the output from net.sf.picard.analysis.CollectAlignmentSummaryMetrics is stored.

###### *.fastq.quality and *fastq.matrix

They contain per-tile quality measures on a space-separated file. Used to create tile mean quality graphs.

###### *.segments

Read length vs proportion of reads. Used to calculate insert size distribution in the summary plots.

<pre>
read_length	proportion_of_reads
0	0.0105
1	0.0007
2	0.0006
3	0.0007
</pre>

##### alignments directory

Contains the results of the alignments *for each sample*. As we see on the listing below, lane 1, barcode id 5 has been aligned in <a href="http://bioinformatics.oxfordjournals.org/content/early/2009/06/08/bioinformatics.btp352.short">SAM</a> and BAM formats. For convenience, to facilitate SNP calling, for instance, a sorted BAM file is also generated.

<pre>
1_110126_B816J0ABXX_5.sam
1_110126_B816J0ABXX_5.bam
1_110126_B816J0ABXX_5-sort.bam
1_110126_B816J0ABXX_5_1_fastq.bam
</pre>

###### *_barcode directories

Those contain fastq files conforming with the naming schema we've seen before. They are the result of the demultiplexing process, where the "unmatched" files contain the reads that have not passed the approximate barcoding matching algorithm:

<pre>
4_110126_B816J0ABXX_1_1.fastq
4_110126_B816J0ABXX_5_2.fastq
4_110126_B816J0ABXX_1_2.fastq
4_110126_B816J0ABXX_6_1.fastq
4_110126_B816J0ABXX_5_1.fastq
4_110126_B816J0ABXX_6_2.fastq

4_110126_B816J0ABXX_unmatched_1.fastq
4_110126_B816J0ABXX_unmatched_2.fastq

4_110126_B816J0ABXX_bc.metrics
SampleSheet-barcodes.cfg
</pre>

*-barcodes.cfg contains a simple mapping between barcode id's and the actual barcode sequence:

<pre>
3 ATCACGA
2 ACTTGAA
9 TAGCTTA
(...)
</pre>

The _bc.metrics file has a plain read distribution for each barcode:

<pre>
2	11594437
3	20247932
9	14390566
unmatched	908420
</pre>

###### How does bcbb barcoding work ?

Barcodes are added to the 3' end of the first sequence. That way, it remains platform-independent and can be easily handled downstream.

On Brad's upcoming <a href="https://github.com/chapmanb/mgh_projects/commit/3387d82f3496025ad13b69e8d9cbb47cf6ee2af9#nglims_paper/nglims_galaxy.tex-P57">paper on ngLIMS</a>, he will explain in some words how demultiplexing works.

The demultiplexing is performed on the cluster by <a href="https://github.com/brainstorm/bcbb/blob/master/nextgen/scripts/barcode_sort_trim.py">barcode_sort_trim.py</a> script.

## FAQ
<pre>
Q: Does it only keep the PF reads?
A: Yes

Q: Does it convert the quality scores to phred-like scores?
A: No, they are encoded in the Illumina 1.3+ format. Looking forward to these being Sanger.

Q: Does if convert dots to Ns?
A: Yes
</pre>

## Requirements

[i1]: http://macnugget.org/projects/publickeys/
[i2]: http://cloudbiolinux.com
[i3]: http://www.biomedcentral.com/1471-2105/11/S12/S4
[i4]: https://github.com/chapmanb/bcbb/blob/master/ec2/biolinux/fabfile.py
[i5]: https://github.com/chapmanb/bcbb/blob/master/ec2/biolinux/data_fabfile.py
[i6]: https://bitbucket.org/afgane/mi-deployment/src
[i7]: http://somethingaboutorange.com/mrl/projects/nose/

### Next gen analysis

* [Picard][3]
* [GATK][4]
* [bowtie][5]
* [samtools][7]
* [snpEff][16]
* [fastx toolkit][6]
* [SolexaQA][6b]
* [matrix2png][6c]

[3]: http://picard.sourceforge.net/
[4]: http://www.broadinstitute.org/gsa/wiki/index.php/The_Genome_Analysis_Toolkit
[5]: http://bowtie-bio.sourceforge.net/
[6]: http://hannonlab.cshl.edu/fastx_toolkit/
[6b]: http://solexaqa.sourceforge.net/
[6c]: http://www.bioinformatics.ubc.ca/matrix2png/
[7]: http://samtools.sourceforge.net/
[16]: http://sourceforge.net/projects/snpeff/

### Processing infrastructure

* RabbitMQ for communication between machines
* LaTeX and pdflatex for report generation
* ps2pdf

### Optional software for generating report graphs

* R with ggplot2, plyr, sqldf libraries.
* [rpy2][11]. Your R needs to be built with shared libraries
  available: `./configure --enable-R-shlib`.

### Python modules installed with the package

* [Biopython][10]
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
