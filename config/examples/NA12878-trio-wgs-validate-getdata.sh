#!/bin/bash
set -eu -o pipefail

# Genome data
wget -c -O NA12878_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
wget -c -O NA12878_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz
wget -c -O NA12891_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194160/ERR194160_1.fastq.gz
wget -c -O NA12891_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194160/ERR194160_2.fastq.gz
wget -c -O NA12892_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194161/ERR194161_1.fastq.gz
wget -c -O NA12892_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194161/ERR194161_2.fastq.gz
