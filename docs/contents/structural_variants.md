# Structural variant calling

## Overview
bcbio can detect larger (>50bp) structural variants like deletions, insertions,
inversions and copy number changes for both germline population and
cancer variant calling

To enable structural variant calling, specify `svcaller` options
in the algorithm section of your configuration

```yaml
- description: Sample
  algorithm:
    svcaller: [lumpy, manta, cnvkit]
```

Split read callers (primary use case - germline WGS sequencing):
- [Lumpy](https://github.com/arq5x/lumpy-sv)
- [Manta](https://github.com/Illumina/manta)
- [WHAM](https://github.com/jewmanchue/wham)
- [DELLY](https://github.com/tobiasrausch/delly)
- we annotate split-read calls with coverage depth using [duphold](https://github.com/brentp/duphold).

Read-depth based CNV callers (primary use case - T/N cancer CNV calling)
- [gatkcnv](https://gatk.broadinstitute.org/hc/en-us/articles/360035531092?id=11682)
- [CNVkit](https://cnvkit.readthedocs.io/en/latest/). 2020-05-13: temporarily off until new release 0.9.7.

## Workflow
This example runs structural variant calling with multiple callers (Lumpy, Manta and CNVkit),
providing a combined output summary file and validation metrics against NA12878 deletions.
It uses the same NA12878 input as the whole genome trio example.

To run the analysis do:
```shell
mkdir -p NA12878-sv-eval
cd NA12878-sv-eval
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-sv-getdata.sh
bash NA12878-sv-getdata.sh
cd work
bcbio_nextgen.py ../config/NA12878-sv.yaml -n 16
```
This is large whole genome analysis and the timing and disk space requirements for the NA12878 trio analysis above apply here as well.

## Paramaters
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

## Validation
- [Validation of germline structural variant detection](https://bcb.io/2014/08/12/validated-whole-genome-structural-variation-detection-using-multiple-callers/) using multiple calling methods to validate against deletions in NA12878. This implements a pipeline that works in tandem with SNP and indel calling to detect larger structural variations like deletions, duplications, inversions and copy number variants (CNVs).
- [Validation of tumor/normal calling](https://bcb.io/2015/03/05/cancerval/) using the synthetic DREAM validation set. This includes validation of additional callers against duplications, insertions and inversions.
- [validation of germline DEL with HuRef benchmark](https://github.com/bcbio/bcbio_validations/blob/master/huref_sv/README.md)

## References
- [Sarwal et al 2020. A comprehensive benchmarking of WGS-based structural variant callers](https://www.biorxiv.org/content/10.1101/2020.04.16.045120v4.full.pdf)
- [Kosugi et al 2019. Comprehensive evaluation of structural variation detection algorithms for whole genome sequencing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1720-5)
- See references to invidivual tools on the citations page.
