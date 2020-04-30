## Configuration

Two configuration files, in easy to write [YAML format](https://en.wikipedia.org/wiki/YAML#Example), specify details about your system and samples to run:
* `bcbio_sample.yaml` Details about a set of samples to process, including input files and analysis options. You configure these for each set of samples to process. This will be the main file prepared for each sample run and the documentation below details techniques to help prepare them.
* `bcbio_system.yaml` High level information about the system, including locations of installed programs like GATK and cores and memory usage (see [Tuning core and memory usage](contents/parallel:tuning%20core%20and%20memory%20usage)). These apply across multiple runs. The automated installer creates a ready to go system configuration file that can be manually edited to match the system. Find the file in the galaxy sub-directory within your installation data location (ie. `/usr/local/share/bcbio-nextgen/galaxy`). To modify system parameters for a specific run, supply [sample or run specific resources](#sample-or-run-specific-resources) in your `bcbio_sample.yaml` file.

Commented [system](https://github.com/bcbio/bcbio-nextgen/blob/master/config/bcbio_system.yaml) and [sample](https://github.com/bcbio/bcbio-nextgen/blob/master/config/bcbio_sample.yaml) example files are available in the `config` directory. The [Example pipelines](contents/intro:example%20pipelines) section contains additional examples of ready to run sample files.

### Automated sample configuration

bcbio-nextgen provides a utility to create configuration files for multiple sample inputs using a base template. Start with one of the [best-practice templates](https://github.com/bcbio/bcbio-nextgen/tree/master/config/templates), or define your own, then apply to multiple samples using the template workflow command:
```shell
bcbio_nextgen.py -w template freebayes-variant project1.csv sample1.bam sample2_1.fq sample2_2.fq
```
* `freebayes-variant` is the name of the standard `freebayes-variant.yaml` input, which the script fetches from GitHub. This argument can also be a path to a locally customized YAML configuration. In both cases, the script replicates the single sample template configuration to all input samples.
* `project1.csv` is a comma separated value file containing sample metadata, descriptions and algorithm tweaks:
    ```
    samplename,description,batch,phenotype,sex,variant_regions
    sample1,ERR256785,batch1,normal,female,/path/to/regions.bed
    sample2,ERR256786,batch1,tumor,,/path/to/regions.bed
    ```
    The first column links the metadata to a specific input file. The template command tries to identify the `samplename` from read group information in a BAM file, or uses the base filename if no read group information is present. For BAM files, this would be the filename without the extension and path (`/path/to/yourfile.bam => yourfile`). For FASTQ files, the template functionality will identify pairs using standard conventions (`_1` and `_2`, including Illumina extensions like `_R1`), so use the base filename without these (`/path/to/yourfile_R1.fastq => yourfile`). Note that paired-end samples sequentially numbered without leading zeros (e.g., `sample_1_1.fastq`, `sample_1_2.fastq`, `sample_2_1.fastq`, `sample_2_2.fastq`, etc., will likely not be parsed correctly; see [issue #1919](https://github.com/bcbio/bcbio-nextgen/issues/1919) for more info). In addition, `.` characters could be problematic, so it's better to avoid this character and use it only as separation for the file extension.

    For [Common Workflow Language (CWL)](cwl) inputs, the first `samplename` column should contain the base filename. For BAM files, this is `your_file.bam`. For fastqs this is `your_file_R1.fastq.gz;your_file_R2.fastq.gz`, separating individual files with a semicolon. By putting paths to the actual locations of the inputs in your `bcbio_system.yaml` input when generating CWL, you can easily move projects between different filesystems.

    The remaining columns can contain:
    * `description` Changes the sample description, originally supplied by the file name or BAM read group, to this value. You can also set the `lane`, although this is less often done as the default sequential numbering works here. 
    * Algorithm parameters specific for this sample. If the column name matches an available [Algorithm parameters](#algorithm-parameters), then this value substitutes into the sample `algorithm`, replacing the defaults from the template. You can also change other information in the BAM read group through the `algorithm` parameters. See [alignment](#alignment) configuration documentation for details on how these map to read group information.
    * metadata key/value pairs. Any columns not falling into the above cases will go into the metadata section. A `ped` specification will allow bcbio to read family, gender and phenotype information from a PED input file and use it for batch, sex and phenotype, respectively. The PED inputs supplement information from the standard template file, so if you specify a value in the template CSV the PED information will no overwrite it. Alternatively, `ped` fields can be specified directly in the metadata as columns. If `family_id` is specified it will be used as the `family_id` for that sample, otherwise `batch` will be used. The `description` column is used as the `individual_id` column and the `phenotype` column will be used for as the `affected` column in the PED format:
        ```
        samplename,description,phenotype,batch,sex,ethnicity,maternal_id,paternal_id,family_id
        NA12878.bam,NA12878,-9,CEPH,female,-9,NA12892,NA12891,NA12878FAM
        ```
    Individual column items can contain booleans (true or false), integers, or lists (separated by semi-colons). These get converted into the expected time in the output YAML file. For instance, to specify a sample that should go into multiple batches:
    ```
    samplename,description,phenotype,batch
    normal.bam,two_normal,normal,Batch1;Batch2
    ```
    For dictionary inputs like [somatic with germline variants](contents/pipelines:somatic%20with%20germline%20variants) setups, you can separate items in a dictionary with colons and double colons, and also use semicolons for lists:
    ```
    samplename,description,phenotype,variantcaller
    tumor.bam,sample1,tumor,germline:freebayes;gatk-haplotype::somatic:vardict;freebayes
    ```
    The name of the metadata file, minus the `.csv` extension, is a short name identifying the current project. The script creates a `project1` directory containing the sample configuration in `project1/config/project1.yaml`.

* The remaining arguments are input BAM or FASTQ files. The script pairs FASTQ files (identified by `_1` and `_2`) and extracts sample names from input BAMs, populating the `files` and `description` field in the final configuration file. Specify the full path to sample files on your current machine.

To make it easier to define your own project specific template, an optional first step is to download and edit a local template. First retrieve a standard template:
```shell
bcbio_nextgen.py -w template freebayes-variant project1
```
This pulls the current GATK best practice variant calling template into your project directory in `project1/config/project1-template.yaml`. Manually edit this file to define your options, then run the full template creation for your samples, pointing to this custom configuration file:
```shell
bcbio_nextgen.py -w template project1/config/project1-template.yaml project1.csv folder/*
```
If your sample folder contains additional BAM or FASTQ files you do not wish to include in the sample YAML configuration, you can restrict the output to only include samples in the metadata CSV with `--only-metadata`. The output will print warnings about samples not present in the metadata file, then leave these out of the final output YAML:
```shell
bcbio_nextgen.py -w template --only-metadata project1/config/project1-template.yaml project1.csv folder/*
```

### Multiple files per sample

In case you have multiple FASTQ or BAM files for each sample you can use `bcbio_prepare_samples.py`. The main parameters are:
* `--out`: the folder where the merged files will be
* `--csv`: the CSV file that is exactly the same as described previously, but having as many duplicate lines for each sample as files to be merged:
```
samplename,description,batch,phenotype,sex,variant_regions
file1.fastq,sample1,batch1,normal,female,/path/to/regions.bed
file2.fastq,sample1,batch1,normal,female,/path/to/regions.bed
file1.fastq,sample2,batch1,tumor,,/path/to/regions.bed
```
An example of usage is:
```shell
bcbio_prepare_samples.py --out merged --csv project1.csv
```
The script will create the `sample1.fastq,sample2.fastq` in the `merged` folder, and a new CSV file in the same folder than the input CSV:`project1-merged.csv`. Later, it can be used for bcbio:
```shell
bcbio_nextgen.py -w template project1/config/project1-template.yaml project1-merged.csv merged/*fastq
```
The new CSV file will look like:
```
samplename,description,batch,phenotype,sex,variant_regions
sample1.fastq,sample1,batch1,normal,female,/path/to/regions.bed
sample2.fastq,sample2,batch1,tumor,,/path/to/regions.bed
```
It supports parallelization the same way `bcbio_nextgen.py` does:
```shell
python $BCBIO_PATH/scripts/utils/bcbio_prepare_samples.py --out merged --csv project1.csv -t ipython -q queue_name -s lsf -n 1
```
See more examples at [parallelize pipeline](parallel).

In case of paired reads, the CSV file should contain all files:
```
samplename,description,batch,phenotype,sex,variant_regions
file1_R1.fastq,sample1,batch1,normal,female,/path/to/regions.bed
file2_R1.fastq,sample1,batch1,normal,female,/path/to/regions.bed
file1_R2.fastq,sample1,batch1,normal,femela,/path/to/regions.bed
file2_R2.fastq,sample1,batch1,normal,female,/path/to/regions.bed
```
The script will try to guess the paired files the same way that `bcbio_nextgen.py -w template` does. It would detect paired files if the difference among two files is only `_R1/_R2` or `-1/-2` or `_1/_2` or `.1/.2`

The output CSV will look like and is compatible with bcbio:
```
samplename,description,batch,phenotype,sex,variant_regions
sample1,sample1,batch1,normal,female,/path/to/regions.bed
```

### Download data from SRA

Use SRA toolkit prefetch/fastq-dump: https://wiki.rc.hms.harvard.edu/display/O2/Aspera+to+download+NCBI+SRA+data

### Sample information

The sample configuration file defines `details` of each sample to process:
```yaml
details:
  - analysis: variant2
    lane: 1
    description: Example1
    files: [in_pair_1.fq, in_pair_2.fq]
    genome_build: hg19
    algorithm:
      platform: illumina
    metadata:
      batch: Batch1
      sex: female
      platform_unit: flowcell-barcode.lane
      library: library_type
```
* `analysis` Analysis method to use [variant2, RNA-seq, smallRNA-seq]
* `lane` A unique number within the project. Corresponds to the `ID` parameter in the BAM read group.
* `description` Unique name for this sample, corresponding to the `SM` parameter in the BAM read group. Required.
* `files` A list of files to process. This currently supports either a single end or two paired-end FASTQ files, or a single BAM file. It does not yet handle merging BAM files or more complicated inputs.
* `genome_build` Genome build to align to, which references a genome keyword in Galaxy to find location build files.
* `algorithm` Parameters to configure algorithm inputs. Options described in more detail below:
  * `platform` Sequencing platform used. Corresponds to the `PL` parameter in BAM read groups. Optional, defaults to `illumina`.
* `metadata` Additional descriptive metadata about the sample:
  * `batch` defines a group that the sample falls in. We perform multi-sample variant calling on all samples with the same batch name. This can also be a list, allowing specification of a single normal sample to pair with multiple tumor samples in paired cancer variant calling (`batch: [MatchWithTumor1, MatchWithTumor2]`).
  * `sex` specifies the sample gender used to correctly prepare X/Y chromosomes. Use `male` and `female` or PED style inputs (1=male, 2=female).
  * `phenotype` stratifies cancer samples into `tumor` and `normal` or case/controls into `affected` and `unaffected`. Also accepts PED style specifications (1=unaffected, 2=affected). CNVkit uses case/control status to determine how to set background samples for CNV calling.
  * `disease` identifies a specific disease name for the sample. Used along with `svprioritize` to help identify gene regions for reporting during analysis with heterogeneity callers like PureCN and TitanCNA. This is primarily for cancer studies and you can narrow genes by disease using inputs like _lung_, _breast_ or _pancreatic_ for different cancer types.
  * `prep_method` A free text description of the method used in sample prep. Used to group together samples during CNV calling for background. This is not required and when not present bcbio assumes all samples in an analysis use the same method.
  * `svclass` defines a classification for a sample for use in SV case/control setups. When set as `control` will put samples into the background samples used for normalization.
  * `ped` provides a [PED phenotype file](http://zzz.bwh.harvard.edu/plink/data.shtml#ped) containing sample phenotype and family information. Template creation uses this to supplement `batch`, `sex` and `phenotype` information provided in the template CSV. GEMINI database creation uses the PED file as input.
  * `platform_unit` -- Unique identifier for sample. Optional, defaults to `lane` if not specified.
  * `library` -- Name of library preparation used. Optional, empty if not present.
  * `validate_batch` -- Specify a batch name to group samples together for preparing validation plots. This is useful if you want to process samples in specific batches, but include multiple batches into the same validation plot.
  * `validate_combine` -- Specify a batch name to combine multiple samples into an additional validation summary. Useful for larger numbers of small samples to evaluate together.

### Upload

The `upload` section of the sample configuration file describes where to put the final output files of the pipeline. At its simplest, you can configure bcbio-nextgen to upload results to a local directory, for example a folder shared amongst collaborators or a Dropbox account. You can also configure it to upload results automatically to a Galaxy instance, to [Amazon S3](https://aws.amazon.com/s3/) or to iRODS. Here is the simplest configuration, uploading to a local directory:
```yaml
upload:
  dir: /local/filesystem/directory
```
General parameters, always required:
* `method` Upload method to employ: `[filesystem, galaxy, s3, irods]`. Defaults to local filesystem.
* `dir` Local filesystem directory to copy to.

Galaxy parameters:
* `galaxy_url` URL of the Galaxy instance to upload to. Upload assumes you are able to access a shared directory also present on the Galaxy machine.
* `galaxy_api_key` User API key to access Galaxy: see the [Galaxy API](https://galaxyproject.org/develop/api/) documentation.
* `galaxy_library` Name of the Galaxy Data Library to upload to. You can specify this globally for a project in `upload` or for individual samples in the sample details section.
* `galaxy_role` Specific Galaxy access roles to assign to the uploaded datasets. This is optional and will default to the access of the parent data library if not supplied. You can specify this globally for a project in `upload` or for individual samples in the sample details section. The [Galaxy Admin](https://galaxyproject.org/data-libraries/#permissions) documentation has more details about roles.

Here is an example configuration for uploading to a Galaxy instance. This assumes you have a shared mounted filesystem that your Galaxy instance can also access:
```yaml
upload:
  method: galaxy
  dir: /path/to/shared/galaxy/filesystem/folder
  galaxy_url: http://url-to-galaxy-instance
  galaxy_api_key: YOURAPIKEY
  galaxy_library: data_library_to_upload_to
```
Your Galaxy `universe_wsgi.ini` configuration needs to have `allow_library_path_paste = True` set to enable uploads.

S3 parameters:
* `bucket` AWS bucket to direct output.
* `folder` A folder path within the AWS bucket to prefix the output.
* `region` AWS region name to use. Defaults to us-east-1
* `reduced_redundancy` Flag to determine if we should store S3 data with reduced redundancy: cheaper but less reliable `[false, true]`

For S3 access credentials, set the standard environmental variables, `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, and `AWS_DEFAULT_REGION` or use [IAM access roles](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/iam-roles-for-amazon-ec2.html) with an instance profile on EC2 to give your instances permission to create temporary S3 access.

iRODS parameters:
* `folder` Full directory name within iRODS to prefix the output.
* `resource` (optional) iRODS resource name, if other than default.

Example configuration:
```yaml
upload:
  method: irods
  dir: ../final
  folder: /irodsZone/your/path/
  resource: yourResourceName
```
Uploads to iRODS depend on a valid installation of the iCommands CLI, and a preconfigured connection through the _iinit_ command.

### Globals

You can define files used multiple times in the `algorithm` section of your configuration in a top level `globals` dictionary. This saves copying and pasting across the configuration and makes it easier to manually adjust the configuration if inputs change:
```yaml
globals:
  my_custom_locations: /path/to/file.bed
details:
  - description: sample1
    algorithm:
      variant_regions: my_custom_locations
  - description: sample2
    algorithm:
      variant_regions: my_custom_locations
```

### Algorithm parameters

The YAML configuration file provides a number of hooks to customize analysis in the sample configuration file. Place these under the `algorithm` keyword.

#### Alignment

* `platform` Sequencing platform used. Corresponds to the `PL` parameter in BAM read groups. Default 'Illumina'.
* `aligner` Aligner to use: [bwa, bowtie, bowtie2, hisat2, minimap2, novoalign, snap, star, tophat2, false] To use pre-aligned BAM files as inputs to the pipeline, set to `false`, which will also skip duplicate marking by default. Using pre-aligned inputs requires proper assignment of BAM read groups and sorting. The `bam_clean` argument can often resolve issues with problematic input BAMs.
* `bam_clean` Clean an input BAM when skipping alignment step. This handles adding read groups, sorting to a reference genome and filtering problem records that cause problems with GATK. Options:
  * `remove_extracontigs` -- Remove non-standard chromosomes (for human, anything that is not chr1-22,X,Y) from the BAM file. This allows compatibility when the BAM reference genome has different contigs from the reference file but consistent ordering for standard chromosomes. Also fixes the read groups in the BAM file as in `fixrg`. This is faster than the full `picard` cleaning option.
  * `fixrg` -- only adjust read groups, assuming everything else in BAM file is compatible.
  * `picard` -- Picard/GATK based cleaning. Includes read group changes, fixing of problematic reads and re-ordering chromosome order to match the reference genome. To fix misencoded input BAMs with non-standard scores, set `quality_format` to `illumina`.
* `bam_sort` Allow sorting of input BAMs when skipping alignment step (`aligner` set to false). Options are coordinate or queryname. For additional processing through standard pipelines requires coordinate sorted inputs. The default is to not do additional sorting and assume pre-sorted BAMs.
* `disambiguate` For mixed or explant samples, provide a list of `genome_build` identifiers to check and remove from alignment. Currently supports cleaning a single organism. For example, with `genome_build: hg19` and `disambiguate: [mm10]`, it will align to hg19 and mm10, run disambiguation and discard reads confidently aligned to mm10 and not hg19. Affects fusion detection when `star` is chosen as the aligner. Aligner must be set to a non false value for this to run.
* `align_split_size`: Increase parallelization of alignment. As of 0.9.8, bcbio will try to determine a useful parameter and you don't need to set this. If you manually set it, bcbio will respect your specification. Set to false to avoid splitting entirely. If set, this defines the number of records to feed into each independent parallel step (for example, 5000000 = 5 million reads per chunk). It converts the original inputs into bgzip grabix indexed FASTQ files, and then retrieves chunks for parallel alignment. Following alignment, it combines all chunks back into the final merged alignment file. This allows parallelization at the cost of additional work of preparing inputs and combining split outputs. The tradeoff makes sense when you have large files and lots of distributed compute. When you have fewer large multicore machines this parameter may not help speed up processing.
* `quality_format` Quality format of FASTQ or BAM inputs [standard, illumina]
* `strandedness` For RNA-seq libraries, if your library is strand specific, set the appropriate flag from [unstranded, firststrand, secondstrand]. Defaults to unstranded. For dUTP marked libraries, firststrand is correct; for Scriptseq prepared libraries, secondstrand is correct.
* `save_diskspace` Remove align prepped bgzip and split BAM files after merging into final BAMs. Helps reduce space on limited filesystems during a run. `tools_off: [upload_alignment]` may also be useful in conjunction with this. `[false, true]`

#### Read trimming

* `trim_reads` Trims low quality or adapter sequences or at the ends of reads using atropos. `adapters` and `custom_trim` specify the sequences to trim. For RNA-seq, it's recommended to leave as False unless running Tophat2. For variant calling, we recommend trimming only in special cases where standard soft-clipping does not resolve false positive problems. Supports trimming with [atropos](https://github.com/jdidion/atropos) or [fastp](https://github.com/OpenGene/fastp). `fastp` is currently not compatible with alignment splitting in variant calling and requires `align_split_size: false`. The old parameter `read_through` defaults to using atropos trimming. `[False, atropos, fastp]`. Default to False.
* `adapters` If trimming adapter read through, trim a set of stock adapter sequences. Allows specification of multiple items in a list, for example [truseq, polya] will trim both TruSeq adapter sequences and polyA tails. polyg trimming removes high quality G stretches present in NovaSeq and NextSeq data. In the small RNA pipeline, bcbio will try to detect the adapter using DNApi. If you set up this parameter, then bcbio will use this value instead. Choices: `[truseq, illumina, nextera, polya, polyx, polyg, nextera2, truseq2]`.
  * nextera2: Illumina NEXTera DNA prep kit from NEB
  * truseq2: SMARTer Universal Low Input RNA Kit
* `custom_trim` A list of sequences to trim from the end of reads, for example: `[AAAATTTT, GGGGCCCC]`
* `min_read_length` Minimum read length to maintain when `read_through` trimming set in `trim_reads`. Defaults to 25.
* `trim_ends` Specify values for trimming at ends of reads, using a fast approach built into fastq preparation. This does not do quality or adapter trimming but will quickly cut off a defined set of values from either the 5' or 3' end of the first and second reads. Expects a list of 4 values: `[5' trim read1, 3' trim read1, 5' trim read2, 3' trim read2]`. Set values to 0 if you don't need trimming (ie. `[6, 0, 12, 0]` will trim 6bp from the start of read 1 and 12bp from the start of read 2. Only implemented for variant calling pipelines.

#### Alignment postprocessing

* `mark_duplicates` Mark duplicated reads [true, false]. If true, will perform streaming duplicate marking with [biobambam's bammarkduplicates or bamsormadup](https://github.com/gt1/biobambam). Uses [samblaster](https://github.com/GregoryFaust/samblaster) as an alternative if you have paired reads and specifying `lumpy` as an `svcaller`. Defaults to true for variant calling and false for RNA-seq and small RNA analyses. Also defaults to false if you're not doing alignment (`aligner: false`).
* `recalibrate` Perform base quality score recalibration on the aligned BAM file, adjusting quality scores to reflect alignments and known variants. Supports both GATK and Sentieon recalibration. Defaults to false, no recalibration. [false, gatk, sentieon]
* `realign` Perform GATK's realignment around indels on the aligned BAM file. Defaults to no realignment since realigning callers like FreeBayes and GATK HaplotypeCaller handle this as part of the calling process. [false, gatk]

#### Coverage information

* `coverage_interval` Regions covered by sequencing. bcbio calculates this automatically from alignment coverage information, so you only need to specify it in the input configuration if you have specific needs or bcbio does not determine coverage correctly. `genome` specifies full genome sequencing, `regional` identifies partial-genome pull down sequencing like exome analyses, and `amplicon` is partial-genome sequencing from PCR amplicon sequencing. This influences GATK options for filtering: we use Variant Quality Score Recalibration when set to `genome`, otherwise we apply cutoff-based soft filters. Also affects copy number calling with CNVkit, structural variant calling and deep panel calling in cancer samples, where we tune regional/amplicon analyses to maximize sensitivity. [genome, regional, amplicon]
* `maxcov_downsample` bcbio downsamples whole genome runs with >10x average coverage to a maximum coverage, avoiding slow runtimes in collapsed repeats and poly-A/T/G/C regions. This parameter specified the multiplier of average coverage to downsample at. For example, _200_ downsamples at 6000x coverage for a 30x whole
genome. Set to _false_ or _0_ to disable downsampling. Current defaults to _false_ pending runtime improvements.
* `coverage_depth_min` Minimum depth of coverage. When calculating regions to call in, bcbio may exclude regions with less than this many reads. It is not a hard filter for variant calling, but rather a guideline for determining callable regions. It's primarily useful when trying to call on very low depth samples. Defaults to 4. Setting lower than 4 will trigger low-depth calling options for GATK.

#### Analysis regions

These BED files define the regions of the genome to analyze and report on. `variant_regions` adjusts regions for small variant (SNP and indel) calling. `sv_regions` defines regions for structural variant calling if different than `variant_regions`. For coverage-based quality control metrics, we first use `coverage` if specified, then `sv_regions` if specified, then `variant_regions`. See the section on [input file preparation](#input-file-preparation) for tips on ensuring chromosome naming in these files match your reference genome. bcbio pre-installs some standard BED files for human analyses. Reference these using the naming schemes described in the [reference data repository](https://github.com/AstraZeneca-NGS/reference_data#capture-region-bed-files).

* `variant_regions` BED file of regions to call variants in.
* `sv_regions` -- A specification of regions to target during structural variant calling. By default, bcbio uses regions specified in `variant_regions` but this allows custom specification for structural variant calling. This can be a pointer to a BED file or special inputs: `exons` for only exon regions, `transcripts` for transcript regions (the min start and max end of exons) or `transcriptsXXXX` for transcripts plus a window of XXXX size around it. The size can be an integer (`transcripts1000`) or exponential (`transcripts1e5`). This applies to CNVkit and heterogeneity analysis.
* `coverage` A BED file of regions to check for coverage and completeness in QC reporting. This can also be a shorthand for a BED file installed by bcbio (see [Structural variant calling](#structural-variant-calling) for options).
* `exclude_regions` List of regions to remove as part of analysis. This allows avoidance of slow and potentially misleading regions. This is a list of the following options:

  * `polyx` Avoid calling variants in regions of single nucleotide stretches greater than 50. These can contribute to long variant calling runtimes when errors in polyX stretches align in high depth to these regions and take a lot of work to resolve. Since we don't expect decent resolution through these types of repeats, this helps avoid extra calculations for assessing the noise. This is an alternative to trimming polyX from the 3' ends for reads with `trim_reads` and `adapters`. Requires an organism with a defined `polyx` file in genome resources. For structural variant calling, adding `polyx` avoids calling small indels for Manta, where these can contribute to long runtimes.
  * `lcr` Avoid calling variants in low complexity regions (LCRs). [Heng Li's variant artifacts paper](https://arxiv.org/abs/1404.0929) provides these regions, which cover ~2% of the genome but contribute to a large fraction of problematic calls due to the difficulty of resolving variants in repetitive regions. Removal can help facilitate comparisons between methods and reduce false positives if you don't need calls in LCRs for your biological analysis. Requires an organism with a defined `lcr` file in genome resources.
  * `highdepth` Remove high depth regions during variant calling, identified by collapsed repeats around centromeres in hg19 and GRCh37 as characterized in the [ENCODE blacklist](https://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/). This is on by default for VarDict and FreeBayes whole genome calling to help with slow runtimes in these regions, and also on for whole genome structural variant calling to avoid false positives from high depth repeats.
  * `altcontigs` Skip calling in unplaced contigs (Un), limit analysis to standard chromosomes -- chr1-22,X,Y,MT for human -- to avoid slowdowns on the additional contigs. By default bcbio calls variants in unplaced but not in alternative contigs (alleles). Alt contig calling is currently not supported.

#### Variant calling

* `variantcaller` Variant calling algorithm. Can be a list of multiple options or false to skip [false, freebayes, gatk-haplotype, haplotyper, platypus, mutect, mutect2, scalpel, tnhaplotyper, tnscope, vardict, varscan, samtools, gatk]
  * Paired (typically somatic, tumor-normal) variant calling is currently supported by vardict, freebayes, mutect2, mutect (see disclaimer below), scalpel (indels only), tnhaplotyper (Sentieon), tnscope (Sentieon) and varscan. See the pipeline documentation on [cancer variant calling](contents/pipelines:cancer%20variant%20calling) for details on pairing tumor and normal samples.
  * You can generate both somatic and germline calls for paired tumor-normal samples using different sets of callers. The pipeline documentation on calling [Somatic with germline variants](contents/pipelines:somatic%20with%20germline%20variants) details how to do this.
  * mutect, a SNP-only caller, can be combined with indels from scalpel or sid. Mutect operates in both tumor-normal and tumor-only modes. In tumor-only mode the indels from scalpel will reflect all indels in the sample, as there is currently no way of separating the germline from somatic indels in tumor-only mode.
* `indelcaller` For the MuTect SNP only variant caller it is possible to add calls from an indelcaller such as scalpel, pindel and somatic indel detector (for Appistry MuTect users only). Currently an experimental option that adds these indel calls to MuTect's SNP-only output. Only one caller supported. Omit to ignore. [scalpel, pindel, sid, false]
* `jointcaller` Joint calling algorithm, combining variants called with the specified `variantcaller`. Can be a list of multiple options but needs to match with appropriate `variantcaller`. Joint calling is only needed for larger input sample sizes (>100
samples), otherwise use standard pooled [population calling](contents/pipelines:population%20calling):
  * `gatk-haplotype-joint` [GATK incremental joint discovery](https://gatkforums.broadinstitute.org/gatk/discussion/3896/the-gatk-reference-model-pipeline-for-incremental-joint-discovery-in-full-detail) with HaplotypeCaller. Takes individual gVCFs called by `gatk-haploype` and perform combined genotyping.
  * `freebayes-joint` Combine freebayes calls using [bcbio.variation.recall](https://github.com/chapmanb/bcbio.variation.recall) with recalling at all positions found in each individual sample. Requires `freebayes` variant calling.
  * `platypus-joint` Combine platypus calls using bcbio.variation.recall with squaring off at all positions found in each individual sample. Requires `platypus` variant calling.
  * `samtools-joint` Combine samtools calls using bcbio.variation.recall with squaring off at all positions found in each individual sample. Requires `samtools` variant calling.
* `joint_group_size` Specify the maximum number of gVCF samples to feed into joint calling. Currently applies to GATK HaplotypeCaller joint calling and defaults to the GATK recommendation of 200. Larger numbers of samples will first get combined prior to genotyping.
* `ploidy` Ploidy of called reads. Defaults to 2 (diploid). You can also tweak specialty ploidy like mitochondrial calling by setting ploidy as a dictionary. The defaults are:
    ```yaml
    ploidy:
      default: 2
      mitochondrial: 1
      female: 2
      male: 1
    ```
* `background` Provide pre-calculated files to use as backgrounds for different processes. Organized as a dictionary with individual keys for different components of the pipeline. You can enter as many or few as needed:
  * `variant` A VCF file with variants to use as a background reference during variant calling. For tumor/normal paired calling use this to supply a panel of normal individuals.
  * `cnv_reference` Background reference file for copy number calling. This can be either a single file for one CNV method or a dictionary for multiple methods. Supports [CNVkit cnn inputs](https://cnvkit.readthedocs.io/en/stable/fileformats.html#copy-number-reference-profile-cnn), [GATK4 HDF5 panel of normals](https://software.broadinstitute.org/gatk/documentation/article?id=11682) and [seq2c](https://github.com/AstraZeneca-NGS/Seq2C) combined mapping plus coverage files:
    ```yaml
    background:
      cnv_reference:
        cnvkit: /path/to/background.cnn
        gatk-cnv: /path/to/background_pon.hdf5
        seq2c: /path/to/background.tsv
    ```

#### Somatic variant calling

* `min_allele_fraction` Minimum allele fraction to detect variants in heterogeneous tumor samples, set as the float or integer __percentage__ to resolve (i.e. 10 = alleles in 10% of the sample). Defaults to 10. Specify this in the tumor sample of a tumor/normal pair. It is percentage, not ratio, it is divided /100.0 when calling vardict!
* `use_lowfreq_filter: false`. When set, forces vardict to report variants with low allelec frequency, useful to call variants in panels with high coverage (>1000x). The default (option is not set to false in the config) is to use low frequency filter, i.e. variants could be underreported (variant VAF is above min_allele_fraction but rejected by the filter).

#### Variant annotation

* `effects` Method used to calculate expected variant effects; defaults to [snpEff](http://snpeff.sourceforge.net/). [Ensembl variant effect predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html) is also available when downloaded using [Customizing data installation](contents/installation:customizing%20data%20installation). [snpeff, vep, false]
* `effects_transcripts` Define the transcripts to use for effect prediction annotation. Options `all`: Standard Ensembl transcript list (the default); `canonical`: Report single canonical transcripts (`-canon` in snpEff, `-pick` in VEP); `canonical_cancer` Canonical transcripts with hand curated changes for more common cancer transcripts (effects snpEff only).
* `vcfanno` Configuration files for [vcfanno](https://github.com/brentp/vcfanno), allowing the application of additional annotations to variant calls. By default, bcbio will try and apply:
  * `gemini` -- External population level annotations from [GEMINI](https://gemini.readthedocs.io). This is only run for human samples with gemini data installed ([Customizing data installation](contents/installation:customizing%20data%20installation)).
  * `somatic` -- Somatic annotations from COSMIC, ClinVar and friends. COSMIC need a custom installation within bcbio ([Customizing data installation](contents/installation:customizing%20data%20installation)). Only added for tumor or tumor/normal somatic calling.
  * `rnaedit` -- RNA editing sites for RNA-seq variant calling runs.
    bcbio installs pre-prepared configuration files in `genomes/build/config/vcfanno` or you can specify the full path to a `/path/your/anns.conf` and optionally an equivalently named `/path/your/anns.lua` file. This value can be a list for multiple inputs.

#### Structural variant calling

* `svcaller` -- List of structural variant callers to use. [lumpy, manta, cnvkit, gatk-cnv, seq2c, purecn, titancna, delly, battenberg]. LUMPY and Manta require paired end reads. cnvkit and gatk-cnv should not be used on the same sample due to incompatible normalization approaches, please pick one or the other for CNV calling.
* `svprioritize` -- Produce a tab separated summary file of structural variants in regions of interest. This complements the full VCF files of structural variant calls to highlight changes in known genes. See the [paper on cancer genome prioritization](https://peerj.com/articles/3166/) for the full details. This can be either the path to a BED file (with `chrom start end gene_name`, see [Input file preparation](#input-file-preparation)) or the name of one of the pre-installed prioritization files:
  * `cancer/civic` (hg19, GRCh37, hg38) -- Known cancer associated genes from [CIViC](https://civic.genome.wustl.edu).
  * `cancer/az300` (hg19, GRCh37, hg38) -- 300 cancer associated genes contributed by [AstraZeneca oncology](https://www.astrazeneca.com/our-focus-areas/oncology.html).
  * `cancer/az-cancer-panel` (hg19, GRCh37, hg38) -- A text file of genes in the AstraZeneca cancer panel. This is only usable for `svprioritize` which can take a list of gene names instead of a BED file.
  * `actionable/ACMG56` -- Medically actionable genes from the [The American College of Medical Genetics and Genomics](http://iobio.io/2016/03/29/acmg56/)
  * `coding/ccds` (hg38) -- [Consensus CDS (CCDS)](https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi) regions with 2bps added to internal introns to capture canonical splice acceptor/donor sites, and multiple transcripts from a single gene merged into a single all inclusive gene entry.
* `fusion_mode` Enable fusion detection in RNA-seq when using STAR (recommended) or Tophat (not recommended) as the aligner. OncoFuse is used to summarise the fusions but currently only supports `hg19` and `GRCh37`. For explant samples `disambiguate` enables disambiguation of `STAR` output [false, true]. This option is deprecated in favor of `fusion_caller`.
* `fusion_caller` Specify a standalone fusion caller for fusion mode. Supports `oncofuse` for STAR/tophat runs, `pizzly` and `ericscript` for all runs. If a standalone caller is specified (i.e. `pizzly` or `ericscript` ), fusion detection will not be performed with aligner. `oncofuse` only supports human genome builds GRCh37 and hg19. `ericscript` supports human genome builds GRCh37, hg19 and hg38 after installing the associated fusion databases ([Customizing data installation](contents/installation:customizing%20data%20installation)).
* `known_fusions` A TAB-delimited file of the format `gene1<tab>gene2`, where `gene1` and `gene2` are identifiers of genes specified under `gene_name` in the attributes part of the GTF file.

#### HLA typing

`hlacaller` -- Perform identification of highly polymorphic HLAs with human build 38 (hg38). The recommended option is `optitype`, using the [OptiType](https://github.com/FRED-2/OptiType) caller. Also supports using the [bwa HLA typing implementation](https://github.com/lh3/bwa/blob/master/README-alt.md#hla-typing) with `bwakit`

#### Validation

bcbio pre-installs standard truth sets for performing validation, and also allows use of custom local files for assessing reliability of your runs:
* `validate` A VCF file of expected variant calls to perform validation and grading of small variants (SNPs and indels) from the pipeline. This provides a mechanism to ensure consistency of calls against a known set of variants, supporting comparisons to
genotyping array data or reference materials.
* `validate_regions` A BED file of regions to evaluate small variant calls in. This defines specific regions covered by the `validate` VCF file.
* `svvalidate` -- Dictionary of call types and pointer to BED file of known regions. For example: `DEL: known_deletions.bed` does deletion based validation of outputs against the BED file.

Each option can be either the path to a local file, or a partial path to a file in the pre-installed truth sets. For instance, to validate an NA12878 run against the [Genome in a Bottle](https://github.com/genome-in-a-bottle) truth set:
```yaml
validate: giab-NA12878/truth_small_variants.vcf.gz
validate_regions: giab-NA12878/truth_regions.bed
svvalidate:
  DEL: giab-NA12878/truth_DEL.bed
```

follow the same naming schemes for small variants, regions and different
structural variant types. bcbio has the following validation materials
for germline validations:
* `giab-NA12878` -- [Genome in a Bottle](https://github.com/genome-in-a-bottle) for NA12878, a Caucasian sample. Truth sets: small_variants, regions, DEL; Builds: GRCh37, hg19, hg38
* `giab-NA24385` -- [Genome in a Bottle](https://github.com/genome-in-a-bottle) for NA24385, an Ashkenazic Jewish sample. Truth sets: small_variants, regions; Builds: GRCh37, hg19, hg38
* `giab-NA24631` -- [Genome in a Bottle](https://github.com/genome-in-a-bottle) for NA24631, a Chinese sample. Truth sets: small_variants, regions; Builds: GRCh37, hg19, hg38
* `giab-NA12878-crossmap` -- [Genome in a Bottle](https://github.com/genome-in-a-bottle) for NA12878 converted to hg38 with CrossMap. Truth sets: small_variants, regions, DEL; Builds: hg38
* `giab-NA12878-remap` -- [Genome in a Bottle](https://github.com/genome-in-a-bottle) for NA12878 converted to hg38 with Remap. Truth sets: small_variants, regions, DEL; Builds: hg38
* `platinum-genome-NA12878` -- [Illumina Platinum Genome](https://www.illumina.com/platinumgenomes/) for NA12878. Truth sets: small_variants, regions; Builds: hg19, hg38

and for cancer validations:
* `giab-NA12878-NA24385-somatic` -- A [sequenced NA12878/NA24385 mixture](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/use_cases/mixtures/UMCUTRECHT_NA12878_NA24385_mixture_10052016/) providing a somatic-like truth set for detecting low frequency events. Build: Truth sets: small_variants, regions. Builds: GRCh37, hg38
* `dream-syn3` -- Synthetic dataset 3 from the [ICGC-TCGA DREAM mutation calling challenge](https://www.synapse.org/#!Synapse:syn312572/wiki/62018). Truth sets: small_variants, regions, DEL, DUP, INV, INS. Builds: GRCh37.
* `dream-syn4` -- Synthetic dataset 4 from the [ICGC-TCGA DREAM mutation calling challenge](https://www.synapse.org/#!Synapse:syn312572/wiki/62018). Truth sets: small_variants, regions, DEL, DUP, INV. Builds: GRCh37.
* `dream-syn3-crossmap` -- Synthetic dataset 3 from the [ICGC-TCGA DREAM mutation calling challenge](https://www.synapse.org/#!Synapse:syn312572/wiki/62018) converted to human build 38 coordinates with CrossMap. Truth sets: small_variants, regions, DEL, DUP, INV, INS. Builds: hg38.
* `dream-syn4-crossmap` -- Synthetic dataset 4 from the [ICGC-TCGA DREAM mutation calling challenge](https://www.synapse.org/#!Synapse:syn312572/wiki/62018) converted to human build 38 coordinates with CrossMap. Truth sets: small_variants, regions, DEL, DUP, INV. Builds: hg38.

For more information on the hg38 truth set preparation see the work on [validation on build 38 and conversion of human build 37 truth sets to build 38](https://bcb.io/2015/09/17/hg38-validation/). The [installation recipes](https://github.com/chapmanb/cloudbiolinux/tree/master/ggd-recipes) contain provenance details about the origins of the installed files.

#### UMIs

Unique molecular identifiers (UMIs) are short random barcodes used to tag individual molecules and avoid amplification biased. Both single cell RNA-seq and variant calling support UMIs. For variant calling, [fgbio](https://github.com/fulcrumgenomics/fgbio) collapses sequencing duplicates for each UMI into a single consensus read prior to running re-alignment and variant calling. This requires `mark_duplicates: true` (the default) since it uses position based duplicates and UMI tags for collapsing duplicate reads into consensus sequences.

To help with preparing fastq files with UMIs bcbio provides a script `bcbio_fastq_umi_prep.py`. This handles two kinds of UMI barcodes:
* Separate UMIs: it converts reads output by an Illumina as 3 files (read 1, read 2, and UMIs).
* Duplex barcodes with tags incorporated at the 5' end of read 1 and read 2

In both cases, these get converted into paired reads with UMIs in the fastq names, allowing specification of `umi_type: fastq_name` in your bcbio YAML configuration. The script runs on a single set of files or autopairs an entire directory of fastq files. To convert a directory with separate UMI files:
```shell
bcbio_fastq_umi_prep.py autopair -c <cores_to_use> <list> <of> <fastq> <files>
```
To convert duplex barcodes present on the ends of read 1 and read 2:
```shell
bcbio_fastq_umi_prep.py autopair -c <cores_to_use> --tag1 5 --tag2 5 <list> <of> <fastq> <files>
```
If you want to prepare your FASTQ files for use with the `umi_type: fastq_name` option and they don't follow the two barcode schemes listed you can pre-transform the FASTQ files yourself by putting `:UMI_yourumisequence` in the read name. Here is an example:
```
@A00574:89:HCLMTDRXX:1:2101:1425:1016:UMI_CGAACGTGTACACG 2:N:0:CTGAAGCT+GGCTCTGA
GTGATATAATTTATTTTCTTAAAATAGCCATGCTGGCTGGAGCCACAGCAGTTTACTCCCAGTTCATTACTCAGCTAACAGACGAAAACCAGT
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
```
For more complex UMI schemes please open up an issue and we can help you write a transformation to prepare your files. We use <https://github.com/vals/umis> to do these more complex transformations.

Configuration options for UMIs:
* `umi_type` The UMI/cellular barcode scheme used for your data. For variant analysis with UMI based consensus calling, supports either `fastq_name` with UMIs in read names or the path to a fastq file with UMIs for each aligned read.
* `correct_umis: [path/to/whitelist_umi.txt]`. For a restricted set of UMIs specify a text file (one UMI per line). UMIs will be corrected with <http://fulcrumgenomics.github.io/fgbio/tools/latest/CorrectUmis.html>

You can adjust the [fgbio default options](https://github.com/bcbio/bcbio-nextgen/blob/8a76c9e546cb79621707082fd763bd643e0e9652/bcbio/ngsalign/postalign.py#L208) by adjusting [resources](#resources). The most common change is modifying the minimum number of reads as input to consensus sequences. This default to 1 to avoid losing reads, 3 is recommended for high depth panels:
```yaml
resources:
  fgbio:
    options: [--min-reads, 3]
```

#### Fast RNA-seq

`transcriptome_fasta` An optional FASTA file of transcriptome sequences to quantitate rather than using bcbio installed transcriptome sequences.

#### smallRNA sequencing

* `adapters` The 3' end adapter that needs to be remove. For NextFlex protocol you can add `adapters: ["4N", "$3PRIME_ADAPTER"]`. For any other options you can use resources: `atropos:options:["-u 4", "-u -4"]`.
* `species` 3 letters code to indicate the species in mirbase classification (i.e. hsa for human).
* `aligner` Currently STAR is the only one tested although bowtie can be used as well.
* `expression_caller` A list of expression callers to turn on: trna, seqcluster, mirdeep2, mirge (read [smallRNA-seq](contents/pipelines:smallrna-seq) to learn how to set up bcbio to run mirge)
* `transcriptome_gtf` An optional GTF file of the transcriptome to for seqcluster.
* `spikein_fasta` A FASTA file of spike in sequences to quantitate.
* `umi_type: 'qiagen_smallRNA_umi'` Support of Qiagen UMI small RNAseq protocol.

#### ChIP/ATAC sequencing

* `peakcaller` bcbio only accepts `[macs2]`
* `aligner` Currently `bowtie2` is the only one tested. `bwa` is also available.
* The `phenotype` and `batch` tags need to be set under `metadata` in the config YAML file. The `phenotype` tag will specify the chip (`phenotype: chip`) and input samples (`phenotype: input`). The `batch` tag will specify the input-chip pairs of samples for example, `batch: pair1`. Same input can be used for different chip samples giving a list of distinct values: `batch: [sample1, sample2]`.
* `chip_method`: currently supporting standard CHIP-seq (TF or broad regions using _chip_) or ATAC-seq (_atac_).
Parameters will change depending on the option to get the best
possible results. Only macs2 supported for now.
* `antibody`: automatically sets peakcalling options tailored for the specific anitbody. Supports _h3f3a_, _h3k27me3_, _h3k36me3_, _h3k4me1_, _h3k79me2_, _h3k9me3_, _h3k9me1_, _h3k9me2_, _h4k20me1_, _h2afz_, _h3ac_, _h3k27ac_, _h3k4me2_, _h3k4me3_, _h3k9ac_, _h3k9me3_.
* `keep_duplicates`: do not remove duplicates before peak calling. Defaults to _False_.
* `keep_multimapped`: do not remove multimappers before peak calling. Defaults to _False_.

You can pass different parameters for `macs2` adding to [resources](#resources):
```yaml
resources:
  macs2:
    options: ["--broad"]
```

#### Methylation

`aligner` supports `bismark`

#### Quality control

* `qc` Allows you to specifically assign quality control modules to run. Generally you want to leave this unset and allow bcbio to run the correct QC metrics for your experiment, or remove specific QC steps you don't want using `tools_off`
([Changing bcbio defaults](#changing-bcbio-defaults)). However, this can allow turning off most of the QC by specifying a single quick running step like `picard`. Available tools are `fastqc`, `samtools`, `coverage`, `picard`, `contamination` (VerifyBamID), `peddy`, `viral`, `damage`, `umi`, `small-rna`, `atropos`, `chipqc`.
* `mixup_check` Detect potential sample mixups. Currently supports [qSignature](https://sourceforge.net/p/adamajava/wiki/qSignature/). `qsignature_full` runs a larger analysis while `qsignature` runs a smaller subset on chromosome 22. [False, qsignature, qsignature_full]
* `kraken` Turn on kraken algorithm to detect possible contamination. You can add `kraken: minikraken` and it will use a minimal database to detect possible [contaminants](https://ccb.jhu.edu/software/kraken/). As well, you can point to a [custom database](https://github.com/DerrickWood/kraken) directory and kraken will use it. You will find the results in the _qc_ directory. You need to use _--datatarget
kraken_ during installation to make the minikraken database available.
* `preseq` Accepts `lc_extrap` or `c_curve`, and runs [Preseq](https://smithlabresearch.org/software/preseq/), a tool that predicts the yield for future experiments. By default, it runs 300 steps of estimation using the segment length of 100000. The default extrapolation limit for `lc_extrap` is 3x of the reads number. You can override the parameters `seg_len`, `steps`, `extrap_fraction` using the [resources](#resources) section:
    ```yaml
    resources:
      preseq:
        extrap_fraction: 5
        steps: 500
        seg_len: 5000
    ```
    And you can also set `extrap` and `step` parameters directly, as well as provide any other command line option via `options`:
    ```yaml
    resources:
      preseq:
        extrap: 10000000
        step: 30000
        options: ["-D"]
    ```
* bcbio uses [MultiQC](https://multiqc.info/) to combine QC output for all samples into a single report file. If you need to tweak configuration settings from bcbio defaults, you can use [resources](#resources). For instance to display read counts with full numbers instead of the default millions:
    ```yaml
    resources:
      multiqc:
        options: ["--cl_config", "'read_count_multiplier: 1'"]
    ```
    or as thousands:
    ```yaml
    resources:
      multiqc:
        options: ["--cl_config", "'{read_count_multiplier: 0.001, read_count_prefix: K}'"]
    ```

#### Post-processing

`archive` Specify targets for long term archival. `cram` removes fastq names and does 8-bin compression of BAM files into [CRAM format](https://www.ebi.ac.uk/ena/software/cram-toolkit). `cram-lossless` generates CRAM files without changes to quality scores or fastq name. Default: [] -- no archiving. Lossy cram has some issues, lossless cram provides pretty good compression relative to BAM, and many machines output binned values now, so `cram-lossless` is what we recommend you use.

#### Changing bcbio defaults

bcbio provides some hints to change default behavior be either turning specific defaults on or off, with `tools_on` and `tools_off`. Both can be lists with multiple options:

* `tools_off` Specify third party tools to skip as part of analysis pipeline. Enables turning off specific components of pipelines if not needed:
  * `gatk4` Use older GATK versions (3.x) for GATK commands like BQSR, HaplotypeCaller and VQSR. By default bcbio includes GATK4 and uses it.
  * `vqsr` turns off variant quality score recalibration for all samples.
  * `bwa-mem` forces use of original `bwa aln` alignment. Without this, we use bwa mem with 70bp or longer reads.
  * `lumpy-genotype` skip genotyping for Lumpy samples, which can be slow in the case of many structural variants.
  * `seqcluster` turns off use of seqcluster tool in srnaseq pipeline.
  * `tumoronly-prioritization` turns off attempted removal of germline variants from tumor only calls using external population data sources like ExAC and 1000 genomes.
  * `vardict_somatic_filter` disables running a post calling filter for VarDict to remove variants found in normal samples. Without `vardict_somatic_filter` in paired analyses no soft filtering of germline variants is performed but all high quality variants pass.
  * `upload_alignment` turns off final upload of large alignment files.
  * `pbgzip` turns off use of bgzip with multiple threads.
  * For quality control, you can turn off any specific tool by adding to `tools_off`. For example, `fastqc` turns off quality control FastQC usage. and `coverage_qc` turns off calculation of coverage statistics with samtools-stats and picard. See the [Methylation](#methylation) docs for details on tools.

* `tools_on` Specify functionality to enable that is off by default:
  * `bcbiornaseq` loads a bcbioRNASeq object for use with [bcbioRNASeq](https://github.com/hbc/bcbioRNASeq).
  * `bnd-genotype` enables genotyping of breakends in Lumpy calls, which improves accuracy but can be slow.
  * `bwa-mem` forces use of bwa mem even for samples with less than 70bp reads.
  * `coverage_perbase` calculates per-base coverage depth for analyzed variant regions.
  * `damage_filter` annotates low frequency somatic calls in INFO/DKFZBias for DNA damage artifacts using [DKFZBiasFilter](https://github.com/eilslabs/DKFZBiasFilter).
  * `gemini` Create a [GEMINI database](https://github.com/arq5x/gemini) of variants for downstream query using the new vcfanno and vcf2db approach.
  * `gemini_allvariants` enables all variants to go into GEMINI, not only those that pass filters.
  * `gemini_orig` Create a [GEMINI database](https://github.com/arq5x/gemini) of variants using the older GEMINI loader. Only works for GRCh37 and hg19.
  * `gvcf` forces gVCF output for callers that support it (GATK HaplotypeCaller, FreeBayes, Platypus). For joint calling using a population of samples, please use _jointcaller_ ([Population calling](contents/pipelines:population%20calling)).
  * `lumpy_usecnv` uses input calls from CNVkit as prior evidence to Lumpy calling.
  * `noalt_calling` call variants only for chr1,,22,X,Y,MT.
  * `qualimap` runs [Qualimap](http://qualimap.bioinfo.cipf.es/) (qualimap uses downsampled files and numbers here are an estimation of 1e7 reads).
  * `qualimap_full` runs Qualimap with full bam files but it may be slow.
  * `svplots` adds additional coverage and summary plots for CNVkit and detected ensemble variants.
  * `tumoronly_germline_filter` applies a `LowPriority` filter to tumor-only calls that match population germline databases. The default is to just apply a tag `EPR` (external prioritization) that flags variants present in external databases. Anything missing a `pass` here is a likely germline.
  * `vcf2db_expand` decompresses and expands the genotype columns in the vcfanno prepared GEMINI databases, enabling standard SQL queries on genotypes and depths.
  * `vqsr` makes GATK try quality score recalibration for variant filtration, even for smaller sample sizes.
  * `vep_splicesite_annotations` enables the use of the MaxEntScan and SpliceRegion plugin for VEP. Both optional plugins add extra splice site annotations.

#### Parallelization

* `nomap_split_size` Unmapped base pair regions required to split analysis into blocks. Creates islands of mapped reads surrounded by unmapped (or N) regions, allowing each mapped region to run in parallel. (default: 250)
* `nomap_split_targets` Number of target intervals to attempt to split processing into. This picks unmapped regions evenly spaced across the genome to process concurrently. Limiting targets prevents a large number of small targets. (default: 200 for standard runs, 20 for CWL runs)

#### Ensemble variant calling

In addition to single method variant calling, we support calling with multiple calling methods and consolidating into a final Ensemble callset.

The recommended method to do this uses a simple majority rule ensemble classifier that builds a final callset based on the intersection of calls. It selects variants represented in at least a specified number of callers:
```yaml
variantcaller: [mutect2, varscan, freebayes, vardict]
ensemble:
  numpass: 2
  use_filtered: false
```
This example selects variants present in 2 out of the 4 callers and does not use filtered calls (the default behavior). Because of the difficulties of producing a unified FORMAT/genotype field across callers, the ensemble outputs contains a mix of outputs from the different callers. It picks a representative sample in the order of specified caller, so in the example above would have a MuTect2 call if present, otherwise a VarScan call if present, otherwise a FreeBayes call. This may require custom normalization scripts during post-processing when using these calls. [bcbio.variation.recall](https://github.com/chapmanb/bcbio.variation.recall) implements this approach, which handles speed and file sorting limitations in the [bcbio.variation](https://github.com/chapmanb/bcbio.variation) approach.

This older approach uses the [bcbio.variation](https://github.com/chapmanb/bcbio.variation) toolkit to perform the consolidation. An example configuration in the
`algorithm` section is:
```yaml
variantcaller: [gatk, freebayes, samtools, gatk-haplotype, varscan]
ensemble:
  format-filters: [DP < 4]
  classifier-params:
    type: svm
  classifiers:
    balance: [AD, FS, Entropy]
    calling: [ReadPosEndDist, PL, PLratio, Entropy, NBQ]
  trusted-pct: 0.65
```
The `ensemble` set of parameters configure how to combine calls from the
multiple methods:
* `format-filters` A set of filters to apply to variants before combining. The example removes all calls with a depth of less than 4.
* `classifier-params` Parameters to configure the machine learning approaches used to consolidate calls. The example defines an SVM classifier.
* `classifiers` Groups of classifiers to use for training and evaluating during machine learning. The example defines two set of criteria for distinguishing reads with allele balance issues and those with low calling support.
* `trusted-pct` Define threshold of variants to include in final callset. In the example, variants called by more than 65% of the approaches (4 or more callers) pass without being requiring SVM filtering.

### Resources

The `resources` section allows customization of locations of programs and memory and compute resources to devote to them:
```yaml
resources:
  bwa:
    cores: 12
    cmd: /an/alternative/path/to/bwa
  samtools:
    cores: 16
    memory: 2G
  gatk:
    jvm_opts: ["-Xms2g", "-Xmx4g"]
  mutect2_filter:
    options: ["--max-events-in-region", "2"]
```
* `cmd` Location of an executable. By default, we assume executables are on the path.
* `cores` Cores to use for multi-proccessor enabled software. This is how many cores will be allocated per job. For example if you are running 10 samples and passed -n 40 to bcbio-nextgen and the step you are running has cores: 8 set, a maximum of five samples will run in parallel, each using 8 cores.
* `jvm_opts` Specific memory usage options for Java software. For memory usage on programs like GATK, specify the maximum usage per core. On multicore machines, that's machine-memory divided by cores. This avoids memory errors when running multiple jobs simultaneously, while the framework will adjust memory up when running multicore jobs.
* `memory` Specify the memory per core used by a process. For programs where memory control is available, like `samtools sort`, this limits memory usage. For other programs this is an estimate of usage, used by [Memory management](contents/parallel:memory%20management) to avoid over-scheduling memory. Always specify this as the memory usage for a single core, and the pipeline handles scaling this when a process uses multiple cores.
* `keyfile` Specify the location of a program specific key file or license server, obtained from a third party software tool. Supports licenses for [novoalign](http://www.novocraft.com/products/novoalign/) and [Sentieon](https://www.sentieon.com/products/). For more complex Sentieon setups this can also be a dictionary of environmental variables:
    ```yaml
    resources:
      sentieon:
        keyfile:
          SENTIEON_LICENSE_SERVER: 100.100.100.100:8888
          SENTIEON_AUTH_MECH: XXX
          SENTIEON_AUTH_DATA: signature
    ```
* `options` Adjust specific command line options for a program. This can be hard to support for many tools due to conflicts with other existing options but is available for some tools: - `mutect2`, `mutect2_filter`: Adjust for Mutect2 calls and filtering.

#### Temporary directory

You also use the resource section to specify system specific parameters like global temporary directories:
```yaml
resources:
  tmp:
    dir: /scratch
```
This is useful on cluster systems with large attached local storage, where you can avoid some shared filesystem IO by writing temporary files to the local disk. When setting this keep in mind that the global temporary disk must have enough space to handle intermediates. The space differs between steps but generally you'd need to have 2x the largest input file per sample and account for samples running simultaneously on multiple core machines.

To handle clusters that specify local scratch space with an environmental variable, bcbio will resolve environmental variables like:
```yaml
resources:
  tmp:
    dir: $YOUR_SCRATCH_LOCATION
```

#### Sample or run specific resources

To override any of the global resource settings in a sample specific manner, you write a resource section within your sample YAML configuration. For example, to create a sample specific temporary directory and pass a command line option to novoalign, write a sample resource specification like:
```yaml
- description: Example
  analysis: variant2
  resources:
    novoalign:
      options: ["-o", "FullNW", "--rOQ"]
    tmp:
      dir: tmp/sampletmpdir
```
To adjust resources for an entire run, you can add this resources specification at the top level of your sample YAML:
```yaml
details:
  - description: Example
resources:
  default:
    cores: 16
```

#### Logging directory

By default, bcbio writes the logging-output directory to `log` in the main directory of the run. You can configure this to a different location in your `bcbio-system.yaml` with:
```yaml
log_dir: /path/to/logs
```

### Input file preparation

Input files for supplementing analysis, like `variant_regions` need to match the specified reference genome. A common cause of confusion is the two chromosome naming schemes for human genome build 37: UCSC-style in hg19 (chr1, chr2) and Ensembl/NCBI style in GRCh37 (1, 2). To help avoid some of this confusion, in build 38 we only support the commonly agreed on chr1, chr2 style.

It's important to ensure that the chromosome naming in your input files match those in the reference genome selected. bcbio will try to detect this and provide helpful errors if you miss it.

To convert chromosome names, you can use [Devon Ryan's collection of chromosome mappings](https://github.com/dpryan79/ChromosomeMappings) as an input to sed. For instance, to convert hg19 chr-style coordinates to GRCh37:
```shell
wget --no-check-certificate -qO- https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh37_UCSC2ensembl.txt \
   | awk '{if($1!=$2) print "s/^"$1"/"$2"/g"}' > remap.sed
sed -f remap.sed original.bed > final.bed
```

### Genome configuration files

Each genome build has an associated `buildname-resources.yaml` configuration file which contains organism specific naming and resource files. bcbio-nextgen expects a resource file present next to the genome FASTA file. [Example genome configuration files](https://github.com/bcbio/bcbio-nextgen/tree/master/config/genomes) are available, and automatically installed for natively supported genomes. Create these by hand to support additional organisms or builds.

The major sections of the file are:
* `aliases` -- Names for third-party programs used as part of the analysis, since naming expectations can differ between software programs.
* `variation` -- Supporting data files for variant analysis. For human analyses, the dbSNP and training files are from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360036212652-Resource-Bundle). These are inputs into the training models for recalibration. The automated [CloudBioLinux](https://github.com/chapmanb/cloudbiolinux) data scripts will download and install these in the variation subdirectory relative to the genome files.
* `rnaseq` -- Supporting data files for RNA-seq analysis. The automated installer and updater handles retrieval and installation of these resources for supported genome builds.
* `srnaseq` -- Supporting data files for smallRNA-seq analysis. Same as in rnaseq, the automated installer and updater handle this for supported genome builds.

By default, we place the `buildname-resources.yaml` files next to the genome FASTA files in the reference directory. For custom setups, you specify an alternative directory in the [resources](#resources) section of your `bcbio_system.yaml` file:
```yaml
resources:
  genome:
    dir: /path/to/resources/files
```

### Reference genome files

The pipeline requires access to reference genomes, including the raw FASTA sequence and pre-built indexes for aligners. The automated installer will install reference files and indexes for commonly used genomes (see the [automated installation](contents/installation:automated) documentation for command line options).

For human genomes, we recommend using build 38 (hg38). This is [fully supported and validated](http://bcb.io/2015/09/17/hg38-validation/) in bcbio, and corrects a lot of issues in the previous build 37. We use the [1000 genomes distribution](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/) which includes HLAs and decoy sequences. For human build 37, GRCh37 and hg19, we use the 1000 genome references provided in the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360036212652-Resource-Bundle). These differ in chromosome naming: hg19 uses chr1, chr2, chr3 style contigs while GRCh37 uses 1, 2, 3. They also differ [slightly in content](https://gatkforums.broadinstitute.org/gatk/discussion/1810/whats-the-difference-between-b37-and-hg19-resources): GRCh37 has masked [Pseudoautosomal regions](https://en.wikipedia.org/wiki/Pseudoautosomal_region) on chromosome Y allowing alignment to these regions on chromosome X.

You can use pre-existing data and reference indexes by pointing bcbio-nextgen at these resources. We use the [Galaxy .loc files](https://galaxyproject.org/admin/data-integration/#set-up-the-loc-file) approach to describing the location of the sequence and index data, as described in [Data requirements](contents/installation:data%20requirements). This does not require a Galaxy installation since the installer sets up a minimal set of `.loc` files. It finds these by identifying the root `galaxy` directory, in which it expects a `tool-data` sub-directory with the `.loc` files. It can do this in two ways:
* Using the directory of your `bcbio-system.yaml`. This is the default mechanism setup by the automated installer and requires no additional work.
* From the path specified by the `galaxy_config` option in your `bcbio-system.yaml`. If you'd like to move your system YAML file, add the full path to your `galaxy` directory here. This is useful if you have a pre-existing Galaxy installation with reference data.

To manually make genomes available to bcbio-nextgen, edit the individual `.loc` files with locations to your reference and index genomes. You need to edit `sam_fa_indices.loc` to point at the FASTA files and then any genome indexes corresponding to aligners you'd like to use (for example: `bwa_index.loc` for bwa and `bowtie2_indices.loc` for bowtie2). The database key names used (like `GRCh37` and `mm10`) should match those used in the `genome_build` of your sample input configuration file.

To remove a reference genome, delete its directory `bcbio/genomes/species/reference` and remove all the records corresponding to that genome from `bcbio/galaxy/tool-data/*.loc` files.

### Adding custom genomes

`bcbio_setup_genome.py` will help you to install a custom genome and apply all changes needed to the configuration files. It needs the genome in FASTA format, and the annotation file in GTF or GFF3 format. It can create index for all aligners used by bcbio. Moreover, it will create the folder _rnaseq_ to allow you run the RNAseq pipeline without further configuration. The `--buildversion` option will write that string to the `version.txt` file, to track from where and which version of a gene build was used.
```shell
bcbio_setup_genome.py -f genome.fa -g annotation.gtf -i bowtie2 star seq -n Celegans -b WBcel135 --buildversion WormBase_34
```
If you want to add smallRNA-seq data files, you will need to add the 3 letters code of mirbase for your genome (i.e hsa for human) and the GTF file for the annotation of smallRNA data. Here you can use the same file than the transcriptome if no other available.
```shell
bcbio_setup_genome.py -f genome.fa -g annotation.gtf -i bowtie2 star seq -n Celegans -b WBcel135 --species cel --srna_gtf another_annotation.gtf --buildversion WormBase_34
```
To use that genome just need to configure your YAML files as:
```yaml
genome_build: WBcel135
```

#### Effects prediction

To perform variant calling and predict effects in a custom genome you'd have to manually download and link this into your installation. First find the snpEff genome build:
```shell
$ snpEff databases | grep Lactobacillus | grep pentosus
Lactobacillus_pentosus_dsm_20314                                Lactobacillus_pentosus_dsm_20314                                              ENSEMBL_BFMPP_32_179            http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_ENSEMBL_BFMPP_32_179.zip
Lactobacillus_pentosus_kca1                                     Lactobacillus_pentosus_kca1                                                   ENSEMBL_BFMPP_32_179            http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_ENSEMBL_BFMPP_32_179.zip
```
then download to the appropriate location:
```shell
$ cd /path/to/bcbio/genomes/Lacto/Lactobacillus_pentosus
$ mkdir snpEff
$ cd snpEff
$ wget http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_ENSEMBL_BFMPP_32_179.zip
$ unzip snpEff_v4_3_ENSEMBL_BFMPP_32_179.zip
$ find . -name "Lactobacillus_pentosus_dsm_20314"
 ./home/pcingola/snpEff/data/Lactobacillus_pentosus_dsm_20314
$ mv ./home/pcingola/snpEff/data/Lactobacillus_pentosus_dsm_20314 .
```
finally add to your genome configuration file (`seq/Lactobacillus_pentosus-resources.yaml`):
```yaml
aliases:
  snpeff: Lactobacillus_pentosus_dsm_20314
```
For adding an organism not present in snpEff, please see this [mailing list discussion](https://groups.google.com/d/msg/biovalidation/LPFBlwVBh5s/AMU7MVvQAwAJ).
