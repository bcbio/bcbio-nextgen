Internals
---------

Overview
~~~~~~~~

.. figure:: images/variant-calling-overview.png
   :align: center
   :alt: Variant calling overview

.. _internals-parallel:

Parallel
~~~~~~~~

bcbio calculates callable regions following alignment using `goleft depth
<https://github.com/brentp/goleft/tree/master/depth>`_. These regions determine
breakpoints for analysis, allowing `parallelization by genomic regions
<http://bcb.io/2013/05/22/scaling-variant-detection-pipelines-for-whole-genome-sequencing-analysis/>`_
during variant calling. Defining non-callable regions allows bcbio to select
breakpoints for parallelization within chromosomes where we won't accidentally
impact small variant detection. The callable regions also supplement the variant
calls to define positions where not called bases are homozygous reference,
as opposed to true no-calls with no evidence. The callable regions plus
variant calls is an alternative to gVCF output which explicitly enumerates
reference calls in the output variant file.

.. figure:: images/parallel-clustertypes.png
   :align: center
   :alt: Parallel approaches

   Overview of cluster types during parallel execution

Somatic tumor only variant calling pipeline with UMIs step by step
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Minor steps (like tabix'ing of vcfs, indexing of bams) and details (full paths)
are omitted.

bcbio.yaml config::

    details:
    - algorithm:
        aligner: bwa
        trim_ends: [2,0,2,0]
        min_allele_fraction: 0.01
        correct_umis: /path/umi.whitelist.txt
        tools_off:
        - gemini
        umi_type: fastq_name
        variantcaller:
        - vardict
        coverage: target.bed
        variant_regions: target.bed
      analysis: variant2
      description: samplex
      files:
      - /path/samplex_1.fq.gz
      - /path/samplex_2.fq.gz
      genome_build: hg38
      metadata:
        phenotype: tumor
    resources:
      fgbio:
        options: [--min-reads, 3]

- step 1: trimming 2p from reads: trim_ends: [2,0,2,0], indexing::

    bgzip --threads 8 -c <(seqtk trimfq -b 2 -e 0 samplex_1.fq.gz) > samplex_1.fq.gz
    bgzip --threads 8 -c <(seqtk trimfq -b 2 -e 0 samplex_2.fq.gz) > samplex_2.fq.gz
    grabix index samplex_2.fq.gz
    grabix index samplex_1.fq.gz

- step 2: alignment with bwa mem, sort and mark duplicates, assign BAM tags (XS=sample, XC=cell, RX=UMI), input: fastq files, output: samplex-sort.bam::

    unset JAVA_HOME && \
    bwa mem -c 250 -M -t 16  -R '@RG\tID:samplex\tPL:illumina\tPU:samplex\tSM:samplex' \
    -v 1 /path/reference/genomes/Hsapiens/hg38/bwa/hg38.fa samplex_1.fq.gz samplex_2.fq.gz  | \
    bamsormadup tmpfile=samplex-sort-sorttmp-markdup inputformat=sam threads=16 \
    outputformat=bam level=0 SO=coordinate | \
    /path/bcbio/anaconda/envs/python2/bin/python /path/bin/umis bamtag - | \
    samtools view -b > samplex-sort.bam

- step 3: correct UMIs with fgbio using the whitelist, input: samplex-sort.bam, output: samplex-sort-umis_corrected.bam::

    unset JAVA_HOME && \
    fgbio -Xms750m -Xmx30g -XX:+UseSerialGC --tmp-dir . --async-io=true --compression=0 \
    CorrectUmis \
    -t XC -m 3 -d 1 -x -U /path/umi.whitelist.txt -i samplex-sort.bam \
    -o samplex-sort-umis_corrected.bam

- step 4: calculate coverage with mosdepth::

    mosdepth -t 16 -F 1804 --no-per-base --by target.bed samplex-rawumi \
    samplex-sort-umis_corrected.bam

- step 5: fgbio GroupReadsByUmi, CallDuplexConsensusReads, FilterConsensusReads with min-reads=3, bam2fastq::

    unset JAVA_HOME && \
    fgbio -Xms750m -Xmx30g -XX:+UseSerialGC --tmp-dir . --async-io=true --compression=0 \
    GroupReadsByUmi \
    --edits=1 --min-map-q=1 -t XC -s paired -i samplex-sort-umis_corrected.bam | \
    fgbio -Xms750m -Xmx30g -XX:+UseSerialGC --tmp-dir . --async-io=true --compression=0 \
    CallDuplexConsensusReads \
    --min-input-base-quality=2 --sort-order=:none: -i /dev/stdin -o /dev/stdout | \
    fgbio -Xms750m -Xmx30g -XX:+UseSerialGC --tmp-dir . --async-io=true --compression=0 \
    FilterConsensusReads --min-reads=3 --min-base-quality=13 --max-base-error-rate=0.1 \
    -r /path/reference/genomes/Hsapiens/hg38/seq/hg38.fa -i /dev/stdin -o /dev/stdout | \
    bamtofastq collate=1 T=samplex-sort-umis_corrected-cumi-1-bamtofastq-tmp \
    F=samplex-sort-umis_corrected-cumi-1.fq.gz F2=samplex-sort-umis_corrected-cumi-2.fq.gz tags=cD,cM,cE gz=1

- step 6: align consensus reads::

    unset JAVA_HOME && bwa mem  -C -c 250 -M -t 16  -R '@RG\tID:samplex\tPL:illumina\tPU:samplex\tSM:samplex' \
    -v 1 /projects/ngs/reference/genomes/Hsapiens/hg38/bwa/hg38.fa \
    samplex-sort-umis_corrected-cumi-1.fq.gz samplex-sort-umis_corrected-cumi-2.fq.gz | \
    samtools sort -@ 16 -m 1G -T samplex-sort-cumi-sorttmp -o samplex-sort-cumi.bam /dev/stdin
    samtools index -@ 16 samplex-sort-cumi.bam samplex-sort-cumi.bam.bai

- step 7: clean variant_regions bed file::

    cat target.bed | grep -v ^track | grep -v ^browser | grep -v ^@ | grep -v ^# | \
    bcbio_python -c 'from bcbio.variation import bedutils; bedutils.remove_bad()' | \
    sort -V -T . -k1,1 -k2,2n > cleaned-target.bed
    cat cleaned-target.bed | bgzip --threads 16 -c > cleaned-target.bed.gz
    tabix -f -p bed cleaned-target.bed.gz
    bedtools merge  -i cleaned-target.bed> cleaned-target-merged.bed
    cat cleaned-target-merged.bed  | bgzip --threads 16 -c > cleaned-target-merged.bed.gz

- step 8: clean coverage bed file (the same in our example)::

    cat target.bed | grep -v ^track | grep -v ^browser | grep -v ^@ | grep -v ^# | \
    iconv -c -f utf-8 -t ascii | sed 's/ //g' | \
    bcbio_python -c 'from bcbio.variation import bedutils; bedutils.remove_bad()' | \
    sort -V -T . -k1,1 -k2,2n > cov-target.bed
    cat cov-target.bed  | bgzip --threads 16 -c > cov-target.bed.gz
    tabix -f -p bed cov-target.bed.gz
    bedtools merge -i cov-target.bed > cov-target-merged.bed
    cat cov-target-merged.bed | bgzip --threads 16 -c > cov-target-merged.bed.gz

- step 9: clean sv regions bed file::

    cat cleaned-target.bed | grep -v ^track | grep -v ^browser | grep -v ^@ | grep -v ^# |  \
    /home/kmhr378/local/bin/bcbio_python -c 'from bcbio.variation import bedutils; bedutils.remove_bad()' | \
    sort -V -T . -k1,1 -k2,2n > \
    svregions-cleaned-target.bed
    cat svregions-cleaned-target.bed | bgzip --threads 16 -c > svregions-cleaned-target.bed.gz

- step10: calculate coverage for 3 bed files with MOSDEPTH::

    export MOSDEPTH_Q0=NO_COVERAGE && export MOSDEPTH_Q1=LOW_COVERAGE && \
    export MOSDEPTH_Q2=CALLABLE && \
    mosdepth -t 16 -F 1804 -Q 1 --no-per-base --by cleaned-target.bed \
    --quantize 0:1:4: samplex-variant_regions samplex-sort-cumi.bam
    mosdepth -t 16 -F 1804  --no-per-base --by svregions-cleaned-target.bed  \
    samplex-sv_regions samplex-sort-cumi.bam
    mosdepth -t 16 -F 1804  --no-per-base --by cov-target.bed samplex-coverage \
    samplex-sort-cumi.bam \
    -T 1,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000

- step11: hts_nim counts::

    hts_nim_tools count-reads -t 16 -F 1804 /path/samplex/counts/fullgenome.bed samplex-sort-cumi.bam > fullgenome-1804-counts.txt
    hts_nim_tools count-reads -t 16 -F 1804 cleaned-target.bed samplex-sort-cumi.bam > cleaned-target-merged-1804-counts.txt

- step12: samtools read statistics::
    samtools stats -@ 16 samplex-sort-cumi.bam > samplex.txt
    samtools idxstats samplex-sort-cumi.bam > samplex-idxstats.txt

- step13: variant calling with vardict (repeated for each alignment chunk)::

    unset R_HOME && unset R_LIBS && unset JAVA_HOME && \
    export VAR_DICT_OPTS='-Xms750m -Xmx3500m -XX:+UseSerialGC -Djava.io.tmpdir=.' && \
    vardict-java -G /path/reference/genomes/Hsapiens/hg38/seq/hg38.fa \
    -N samplex -b samplex-sort-cumi.bam -c 1 -S 2 -E 3 -g 4 --nosv --deldupvar -Q 10 -F 0x700 \
    samplex-chr5_0_x-unmerged-regions-regionlimit.bed -f 0.0025 | \
    teststrandbias.R | \
    var2vcf_valid.pl -A -N samplex -E -f 0.0025 | grep -v ^##contig | \
    bcftools annotate -h samplex-chr5_0_x-contig_header.txt | \
    bcftools filter -i 'QUAL >= 0' | \
    bcftools filter --soft-filter 'LowFreqBias' --mode '+' -e  'FORMAT/AF[0:0] < 0.02 && \
    FORMAT/VD[0] < 30 && INFO/SBF < 0.1 && INFO/NM >= 2.0' | \
    awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", $4) } {print}' | \
    awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", $5) } {print}' |  \
    awk -F$'\t' -v OFS='\t' '$1!~/^#/ && $4 == $5 {next} {print}' | \
    vcfstreamsort | bgzip -c > samplex-chr5_0_x.vcf.gz
    zgrep ^# samplex-chr5_0_x.vcf.gz > samplex-chr5_0_x-fixheader-header.vcf
    unset JAVA_HOME && \
    picard FixVcfHeader HEADER=samplex-chr5_0_x-fixheader-header.vcf \
    INPUT=samplex-chr5_0_x.vcf.gz \
    OUTPUT=samplex-chr5_0_x-fixheader.vcf.gz

- step14: gather vcfs::

    unset JAVA_HOME && \
    gatk --java-options '-Xms681m -Xmx3181m -XX:+UseSerialGC -Djava.io.tmpdir=.\
    GatherVcfs -I samplex-files.list -O samplex.vcf.gz

- step15: annotate with snpEff::

    unset JAVA_HOME && \
    snpEff -Xms750m -Xmx29g -Djava.io.tmpdir=. eff \
    -dataDir /path/reference/genomes/Hsapiens/hg38/snpeff \
    -hgvs -cancer -noLog -i vcf -o vcf -csvStats samplex-effects-stats.csv \
    -s samplex-effects-stats.html GRCh38.86 samplex.vcf.gz | \
    bgzip --threads 16 -c > samplex-effects.vcf.gz

- step16: annotate with vcfanno::

    vcfanno -p 16 dbsnp.conf samplex-effects.vcf.gz | \
    bcftools reheader -h samplex-effects-annotated-sample_header.txt | \
    bcftools view | bgzip -c > samplex-effects-annotated.vcf.gz
    tabix -f -p vcf samplex-effects-annotated.vcf.gz
    vcfanno -p 16 samplex-effects-annotated-annotated-somatic-combine.conf \
    samplex-effects-annotated.vcf.gz | bgzip -c > samplex-effects-annotated-annotated-somatic.vcf.gz
    tabix -f -p vcf samplex-effects-annotated-annotated-somatic.vcf.gz
    cat samplex-effects-annotated-annotated-somatic-priority.tsv  | bgzip --threads 16 -c > \
    samplex-effects-annotated-annotated-somatic-priority.tsv.gz
    tabix -f -0 -c '#' -s 1 -b 2 -e 3 samplex-effects-annotated-annotated-somatic-priority.tsv.gz

scRNA-seq: Single cell RNA-seq barcode counting pipeline step by step
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Minor steps and exact file locations are omitted.
* = repeated for every sample in sample_barcodes.csv

bcbio.yaml config::

    details:
    - algorithm:
        cellular_barcode_correction: 1
        minimum_barcode_depth: 1000
        sample_barcodes: /path/project/config/barcodes.csv
        transcriptome_fasta: /path/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.fa
        transcriptome_gtf: /path/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.gtf
        umi_type: harvard-indrop-v3
      analysis: scRNA-seq
      description: project
      files:
      - /path/project_1.fq.gz
      - /path/project_2.fq.gz
      - /path/project_3.fq.gz
      - /path/project_4.fq.gz
      genome_build: mm10
    fc_name: sc-mouse
    upload:
      dir: /path/project/final

- step1: parse barcode information from reads 2,3,4 to the fastq read name,
  CELL_[value]:UMI_[value]:sample_[value]. umis supports many protocols, however,
  the downside is speed - this step can take up to 3-4 days::

    python umis fastqtransform \
    --separate_cb umis/harvard-indrop-v3-transform.json --cores 16  \
    project_1.fq.gz project_2.fq.gz project_3.fq.gz project_4.fq.gz | \
    seqtk seq -L 20 - | gzip > project.umitransformed.fq.gz

- step2: create a fastq file for each sample::

    python umis demultiplex_samples \
    --nedit 1 --barcodes sample_barcodes.csv \
    --out_dir demultiplexed project.umitransformed.fq.gz

- step3*::

    python umis cb_filter \
    --cores 16 --bc1 harvard-indrop-v3-cb1.txt.gz --nedit 1 \
    --bc2 harvard-indrop-v3-cb2.txt.gz demultiplexed/[sample-barcodeAATTTTT].fq  | \
    gzip -c > project-sample_barcode.filtered.fq.gz

- step 4*: create cellular barcode histogram, also creates cb-histogram-filtered.txt for cells
  with nreads > minimum_barcode_depth::

    python umis cb_histogram project-sample_barcode.filtered.fq.gz > cb-histogram.txt

- step 5: create index genome for rapmap::

    rapmap quasiindex -k 31 -i mm10 -t mm10/ref-transcripts.fa

- step 6*. align reads with rapmap::

    rapmap quasimap -t 16 -i mm10 \
    -r <(gzip -cd project-sample_barcode.filtered.fq.gz) | \
    samtools sort -@ 16 -m 1G  -T project-sample_barcode-sorttmp \
    -o project-sample_barcode.bam /dev/stdin
    samtools index -@ 16 project-sample_barcode.bam project-sample_barcode.bam.bai

- step 7*: count transcripts::

    python umis fasttagcount --cb_cutoff 1000 \
    --genemap ref-transcripts-tx2gene.tsv
    --cb_histogram project-sample_barcode/cb-histogram.txt \
    --umi_matrix project-sample_barcode-dupes.mtx.full \
    project-sample_barcode.bam project-sample_barcode.mtx.full

    python umis sparse project-sample_barcode.mtx.full \
    project-sample_barcode.mtx

    python umis sparse project-sample_barcode-dupes.mtx.full \
    project-sample_barcode-dupes.mtx

- step 8. Concatenate all cb-histogram-filtered.txt files::

    cat project-[all-barcodes]/cb-histogram-filtered.txt > cb-histogram.txt

Tests
~~~~~
To run bcbio automated tests, install bcbio and clone bcbio master repository. You are testing your installation with tests provided in bcbio-nextgen/tests::

    which bcbio_nextgen.py
    cd bcbio-nextgen/tests
    ./run_tests.sh > tests.out

Tests are in integration/*.py. Each test has a set or marks. Marks are listed in pytest.ini. The mark defines how many tests to select. By default (just running plain ./run_tests.sh), it is speed1 = 11 tests.

Profiling
~~~~~~~~~
Profiling (tracking CPU, memory, IO usage) could help to optimize resource
usage of bcbio, especially when running  on a server or AWS instance.
Sometimes running a bcbio project with 32 cores is just 10% more efficient
than with 16 cores, because a particular configuration might have memory or
IO related bottlenecks.

- step 1. install and start sysstat deamon: http://www.leonardoborda.com/blog/how-to-configure-sysstatsar-on-ubuntudebian/.
- step 2. Create a cron job to gather system statistics every minute or two.
- step 3. Before bcbio start, drop system memory caches. Otherwise memory usage statistic might be misleading::

    sudo su
    echo 1 > /proc/sys/vm/drop_caches

- step 4. Record bcbio project start and stop time (`date`)
- step 5. collect usage statistics::

    # CPU load
    sar -q -s $start -e $end |  awk '{print $1","$4}'  | sed 1d | sed 1d > cpu.csv
    # memory
    sar -r -s $start -e $end | awk '{print $5}' | sed 1d | sed 1d > mem.csv
    # IO
    sar -b -s $start -e $end | awk '{print $5","$6}' | sed 1d | sed 1d > io.csv
    paste -d "," cpu.csv mem.csv io.csv > usage.csv

    # Example of usage.csv, man sar for more options
    head usage.csv

    23:07:01,ldavg-1,%memused,bread/s,bwrtn/s
    23:08:01,1.77,3.10,23238.66,20204.10
    23:09:01,4.34,24.28,208650.45,26270.58
    23:10:01,11.09,25.67,0.13,15.46
    23:11:01,13.56,27.00,4.27,21.99
    23:12:01,15.44,29.63,26.52,2749.22
    23:13:01,15.16,29.75,42.93,27.06
    23:14:01,16.94,30.54,205.26,2740.95
    23:15:01,15.76,30.57,28.92,2751.62
    23:16:01,15.77,30.88,6.13,33.59

- step 6. Overlap profiling results with bcbio-nextgen-commands.log to investigate the
    performance of particular steps.
