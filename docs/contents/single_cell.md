# Counting cells with bcbio for inDrops3 data: proto-SOP

Bcbio installation paths in this workflow correspond to [O2 bcbio installation](https://wiki.rc.hms.harvard.edu/display/O2).
Adjust to bcbio installation you are working with.

## 1. Check reference genome and transcriptome - is it a mouse project?
- mm10 reference genome: /n/shared_db/bcbio/biodata/genomes/Mmusculus/mm10
- transcriptome_fasta: /n/shared_db/bcbio/biodata/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.fa
- transcriptome_gtf: /n/shared_db/bcbio/biodata/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.gtf

## 2. Create bcbio project structure in /scratch
```
mkdir sc_mouse
cd sc_mouse
mkdir config input final work
```

## 3. Prepare fastq input in sc_mouse/input
- some FC come in 1..4 lanes, merge lanes for every read:
```
cat lane1_r1.fq.gz lane2_r1.fq.gz > project_1.fq.gz
cat lane1_r2.fq.gz lane2_r2.fq.gz > project_2.fq.gz
```
- cat'ing gzip files sounds ridiculous, but works for the most part, for purists:
```
zcat KM_lane1_R1.fastq KM_lane2_R1.fastq.gz | gzip > KM_1.fq.gz
```

- some cores send bz2 files not gz
```
bunzip2 *.bz2
cat *R1.fastq | gzip > sample_1.fq.gz
```

- some cores produce R1,R2,R3,R4, others R1,R2,I1,I2, rename them
```
bcbio_R1 = R1 = 86 or 64 bp transcript read
bcbio_R2 = I1 = 8 bp part 1 of cell barcode
bcbio_R3 = I2 = 8 bp sample (library) barcode
bcbio_R4 = R2 = 14 bp = 8 bp part 2 of cell barcode + 6 bp of transcript UMI
```
- files in sc_mouse/input should be (KM here is project name):
```
KM_1.fq.gz
KM_2.fq.gz
KM_3.fq.gz
KM_4.fq.gz
```

## 4. Specify sample barcodes
Sample barcodes should be in *sc_mouse/config/sample_barcodes.csv*.    

Check out if the sample barcodes provided match the actual barcodes in the data.

```shell
gunzip -c FC_X_3.fq.gz | awk '{if(NR%4 == 2) print $0}' | head -n 400000 | sort | uniq -c | sort -k1,1rn | awk '{print $2","$1}' | head

AGGCTTAG,112303
ATTAGACG,95212
TACTCCTT,94906
CGGAGAGA,62461
CGGAGATA,1116
CGGATAGA,944
GGGGGGGG,852
ATTAGACC,848
ATTAGCCG,840
ATTATACG,699
```

Sometimes you need to reverse complement sample barcodes:
```
cat barcodes_original.csv | awk -F ',' '{print $1}' | tr ACGTacgt TGCAtgca | rev
```

sample_barcodes.csv
```
TCTCTCCG,S01
GCGTAAGA,S02
CCTAGAGT,S03
TCGACTAG,S04
TTCTAGAG,S05
```

## 5. Create bcbio yaml config file
*sc_mouse/config/sc-mouse.yaml*:

```
details:
- algorithm:
    cellular_barcode_correction: 1
    minimum_barcode_depth: 1000
    sample_barcodes: /full/path/sc_mouse/config/sample_barcodes.csv
    transcriptome_fasta: /n/shared_db/bcbio/biodata/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.fa
    transcriptome_gtf: /n/shared_db/bcbio/biodata/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.gtf
    umi_type: harvard-indrop-v3
  analysis: scRNA-seq
  description: PI_name
  files:
  - /full/path/sc_mouse/input/KM_1.fq.gz
  - /full/path/sc_mouse/input/KM_2.fq.gz
  - /full/path/sc_mouse/input/KM_3.fq.gz
  - /full/path/sc_mouse/input/KM_4.fq.gz
  genome_build: mm10
  metadata: {}
fc_name: sc-mouse
upload:
  dir: /full/path/sc_mouse/final
```
Use `cd sc_mouse/input; readlink -f *` to grab full path to each file and paste into yaml.

## 6. Create a batch script
*sc_mouse/config/bcbio.sh*:

```shell
#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=priority        # Partition (queue)
#SBATCH --time=10-00:00             # Runtime in D-HH:MM format
#SBATCH --job-name=km            # Job name
#SBATCH -c 20
#SBATCH --mem-per-cpu=5G            # Memory needed per CPU
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

bcbio_nextgen.py ../config/sc-mouse.yaml -n 20
```
- most projects take < 5days, but some large 4 lane could take more, like 7-8

## 7. Run bcbio

```shell
cd sc_mouse_work
sbatch ../config/bcbio.sh
```

## 1a. (Optional).
If you care, download fresh transcriptome annotation from Gencode (https://www.gencodegenes.org/mouse/)
(it has chrom names with chr matching mm10 assembly).
```
cd sc_mouse/input
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz
gunzip gencode.vM23.annotation.gtf.gz
gffread -g /n/shared_db/bcbio/biodata/genomes/Mmusculus/mm10/seq/mm10.fa gencode.vM23.annotation.gtf -x gencode.vM23.annotation.cds.fa
```
update sc_mouse/config/sc_mouse.yaml:
```
transcriptome_fasta: gencode.vM23.annotation.cds.fa
transcriptome_gtf: gencode.vM23.annotation.gtf
```
## References
- [Indrops3 library structure](https://singlecellcore.hms.harvard.edu/resources)
- [Even shorter guide](https://github.com/bcbio/bcbio-nextgen/blob/master/config/templates/indrop-singlecell.yaml)
- [Much more comprehensive guide](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/01_bcbio_run.md)
