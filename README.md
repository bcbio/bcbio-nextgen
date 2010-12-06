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

## Code structure

Two main scripts drive the automation of the process:

* `scripts/illumina_finished_msg.py` -- Sits on a machine where sequencing
  runs are dumped. It checks for new results, reporting to a RabbitMQ messaging
  queue whenever a new run is finished.
* `scripts/analyze_finished_sqn.py` -- Continuously running server script on
  the Galaxy analysis machine. When new results are reported in the messaging queue,
  this copies over the relevant files and kicks off an automated analysis.

The scripts involved in the actual processing:

* `scripts/automated_initial_analysis.py` -- Drives the high level analysis of
  sequencing lanes based on information specified through the Galaxy LIMS system
* `scripts/upload_to_galaxy.py` -- Handles storing and uploading Fastq,
  alignment, analysis and summary files to Galaxy.
* `scripts/align_summary_report.py` -- Produces a PDF summary file with
  statistics on alignments, duplicates, GC distribution, quality scores,
  and other metrics of interest.
 
System specific information is specified in YAML configuration files:

* `config/post_process.yaml` -- The main configuration file containing Galaxy details,
  program commandlines and customization for processing algorithms.
* `config/transfer_info.yaml` -- Configuration on the sequencing machine, specifying where
  to check for new results.

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

[10]: http://biopython.org
[11]: http://rpy.sourceforge.net/rpy2.html
[12]: http://code.google.com/p/pysam/
[13]: http://www.makotemplates.org/
[14]: http://pyyaml.org/
[15]: http://code.google.com/p/py-amqplib
