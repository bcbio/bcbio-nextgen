# Citations

If you use bcbio in your work, please cite it via Zenodo to help us track how and where it is being used:
<https://doi.org/10.5281/zenodo.3564938>

Please cite individual tools in the manuscript methods, not just bcbio.
To promote citation, we are maitaining easy-to-paste reference lists below. 
Feel free to contribute, if you see a tool that is used through bioconda, but is not included here.

Try to look at the particular tool's documentation, article, and github issues before submitting
an issue to bcbio. This might help to discriminate early between bcbio issues and tools' issues.
The authors of the tool might be of better help to fix the issue. Raising the issue in tool's
github and referencing it in bcbio github issue could speed up its resolution.

## Variant calling
Overall, the parameters of our variant calling workflows are based on GATK best practices (https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows), contributions from bcbio community (https://github.com/bcbio/bcbio-nextgen) and our validations (https://github.com/bcbio/bcbio_validations/).

## Read alignment
We align reads with `bwa mem` [1], using samtools [2], and sambamba [3], to sort bam files and mark duplicate reads.

1. Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. 2013 arXiv:1303.3997. https://github.com/lh3/bwa.
2. Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, and 1000 Genome Project Data Processing Subgroup, The Sequence alignment/map (SAM) format and SAMtools, Bioinformatics (2009) 25(16) 2078-9 [19505943]. https://github.com/samtools/.
3. A. Tarasov, A. J. Vilella, E. Cuppen, I. J. Nijman, and P. Prins. Sambamba: fast processing of NGS alignment formats. Bioinformatics, 2015. https://github.com/biod/sambamba.

## Interval arithmetics.
We use bedtools[1] and pybedtools[2] to wok with genomic intervals

1. Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26, 6, pp. 841–842.
2. Dale RK, Pedersen BS, and Quinlan AR. Pybedtools: a flexible Python library for manipulating genomic datasets and annotations. Bioinformatics (2011). doi:10.1093/bioinformatics/btr539

## Quality control
We run many tools to gather QC metrics:
- peddy [1]
- verifybamid (https://genome.sph.umich.edu/wiki/VerifyBamID)
- DKFZ bias filter(https://github.com/DKFZ-ODCF/DKFZBiasFilter)
- fastqc (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- qualimap (http://qualimap.bioinfo.cipf.es/)
- samtools (https://github.com/samtools/)
- bcftools (http://www.htslib.org/doc/bcftools.html)

We aggregate all metrics in a single QC report with multiqc [2].

1. Pedersen BS, Quinlan AR. Who's Who? Detecting and Resolving Sample Anomalies in Human DNA Sequencing Studies with Peddy. Am J Hum Genet. 2017;100(3):406‐413. doi:10.1016/j.ajhg.2017.01.017 https://github.com/brentp/peddy
2. Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PMID: 27312411; PMCID: PMC5039924.https://multiqc.info/

## Coverage and callable regions
We calculate coverage using mosdepth [1] and calculate callable regions based on real coverage and bed files provided.

1. Pedersen BS, Quinlan AR. Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics. 2018;34(5):867‐868. doi:10.1093/bioinformatics/btx699. https://github.com/brentp/mosdepth.

## SNP and indels in germline (WES, WGS, gene panels)
We support variant calling in germline with:
- gatk4x (https://github.com/broadinstitute/gatk/)
- gatk3.8x (https://console.cloud.google.com/storage/browser/gatk-software/package-archive)
- freebayes (https://github.com/ekg/freebayes)
- samtools (https://github.com/samtools/)
- octopus (https://github.com/luntergroup/octopus).

We support:
- single sample variant calling
- batch variant calling
- population variant calling (joint genotyping).

## Structural and copy number variants in germline (WGS data)
We call structural variants with 
- manta [1]
- lumpy [2]
- delly [3]
- wham [4]

We annotate structural variant calls with coverage information using duphold (https://github.com/brentp/duphold)

1. Chen, X. et al. (2016) Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. Bioinformatics, 32, 1220-1222. doi:10.1093/bioinformatics/btv710. https://github.com/Illumina/manta
2. Layer RM, Chiang C, Quinlan AR, Hall IM. LUMPY: a probabilistic framework for structural variant discovery. Genome Biol. 2014;15(6):R84. Published 2014 Jun 26. doi:10.1186/gb-2014-15-6-r84. https://github.com/arq5x/lumpy-sv
3. Rausch T, Zichner T, Schlattl A, Stütz AM, Benes V, Korbel JO. DELLY: structural variant discovery by integrated paired-end and split-read analysis. Bioinformatics. 2012;28(18):i333‐i339. doi:10.1093/bioinformatics/bts378. https://github.com/dellytools/delly
4. Kronenberg ZN, Osborne EJ, Cone KR, et al. Wham: Identifying Structural Variants of Biological Consequence. PLoS Comput Biol. 2015;11(12):e1004572. Published 2015 Dec 1. doi:10.1371/journal.pcbi.1004572. https://github.com/zeeev/wham

## Somatic small variants
We call somatic variants in tumor only or tumor/normal mode with:
- mutect2 (https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)
- vardict (https://github.com/AstraZeneca-NGS/VarDict)
- strelka2 (https://github.com/Illumina/strelka)
- varscan2 (Koboldt DC, Zhang Q, Larson DE, Shen D, McLellan MD, Lin L, Miller CA, Mardis ER, Ding L, & Wilson RK (2012). VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing. Genome Research PMID: 22300766)

We support ensemble approach to combine somatic calls from several callers (https://github.com/bcbio/bcbio.variation.recall)

## Somatic copy number variants
We use 
- gatk-cnv (https://gatkforums.broadinstitute.org/gatk/discussion/9143/how-to-call-somatic-copy-number-variants-using-gatk4-cnv)
- pureCN (https://bioconductor.org/packages/release/bioc/html/PureCN.html)
- seq2c (https://github.com/AstraZeneca-NGS/Seq2C)
- titanCNA (https://github.com/gavinha/TitanCNA)
- cnvkit

## Variant annotation
We annotate variants with
- VEP [1]
- snpEff [2]
- using vcfanno [3] (https://github.com/brentp/vcfanno) and many annotation sources:
  - gnomad (https://gnomad.broadinstitute.org/)
  - topmed (https://bravo.sph.umich.edu/freeze5/hg38/)
  - cosmic (https://cancer.sanger.ac.uk/cosmic)
  - dbsnp  (https://www.ncbi.nlm.nih.gov/snp/)
  - dbnsfp (https://sites.google.com/site/jpopgen/dbNSFP)
  - clinvar (https://www.ncbi.nlm.nih.gov/clinvar/)

We create a gemini database [4] as output (https://gemini.readthedocs.io/en/latest/). We support any internal vcf or bed based annotation (internal frequency database) via vcfanno.

1. McLaren W, Gil L, Hunt SE, et al. The Ensembl Variant Effect Predictor. Genome Biol. 2016;17(1):122. Published 2016 Jun 6. doi:10.1186/s13059-016-0974-4. https://useast.ensembl.org/info/docs/tools/vep/index.html
2. A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.", Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (Austin). 2012 Apr-Jun;6(2):80-92. PMID: 22728672].
3. Pedersen BS, Layer RM, Quinlan AR. Vcfanno: fast, flexible annotation of genetic variants. Genome Biol. 2016;17(1):118. Published 2016 Jun 1. doi:10.1186/s13059-016-0973-5.
4. Paila U, Chapman BA, Kirchner R, Quinlan AR. GEMINI: integrative exploration of genetic variation and genome annotations. PLoS Comput Biol. 2013;9(7):e1003153. doi:10.1371/journal.pcbi.1003153.

## bulk RNA-seq
We align with STAR[1] or hisat[2], quantify transcripts with salmon[3], kallisto[4], stringtie [5].
1. Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15‐21. doi:10.1093/bioinformatics/bts635
2. Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol. 2019;37(8):907‐915. doi:10.1038/s41587-019-0201-4
3. Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware quantification of transcript expression. Nat Methods. 2017;14(4):417‐419. doi:10.1038/nmeth.4197
4. Bray NL, Pimentel H, Melsted P, Pachter L. Near-optimal probabilistic RNA-seq quantification [published correction appears in Nat Biotechnol. 2016 Aug 9;34(8):888]. Nat Biotechnol. 2016;34(5):525‐527. doi:10.1038/nbt.3519
5. Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT, Salzberg SL. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nat Biotechnol. 2015;33(3):290‐295. doi:10.1038/nbt.3122 

## Fusion calling - RNA-seq
We call fusions with oncofuse[1], pizzly[2], ericscript[3], arriba[4]. 
4. Uhrig S. Arriba - Fast and accurate gene fusion detection from RNA-Seq data 2019. Available from: https://github.com/suhrig/arriba. 

## ATAC-seq
We are using ataqv[1] for quality control.

1. Orchard P, Kyono Y, Hensley J, Kitzman JO, Parker SCJ. Quantification, Dynamic Visualization, and Validation of Bias in ATAC-Seq Data with ataqv. Cell Syst. 2020;10(3):298‐306.e4. doi:10.1016/j.cels.2020.02.009 https://www.cell.com/cell-systems/pdfExtended/S2405-4712(20)30079-X https://github.com/ParkerLab/ataqv 

## small RNA-seq

Data was analyzed with bcbio-nextgen (<https://github.com/bcbio/bcbio-nextgen>) using piDNA to detect the adapter, cutadapt to remove it, STAR/bowtie to align against the genome and seqcluster to detect small RNA transcripts. miRNAs were detected using miraligner tool with miRBase as the reference miRNA database. tRNA profiles were detected using tdrmapper tool. mirdeep2 was used for discovery of novel miRNAs. FastQC was used for QC metrics and multiqc for reporting.

Download BIB format: <https://github.com/bcbio/bcbio-nextgen/tree/master/docs/contents/misc/bcbio-smallrna.bib>

### Tools

* Tsuji J, Weng Z. (2016) DNApi: A De Novo Adapter Prediction Algorithm for Small RNA Sequencing Data. 11(10):e0164228. <https://doi.org/10.1371/journal.pone.0164228>
* Andrews, S. (2010). FastQC: A quality control tool for high throughput sequence data. Bioinformatics. doi:citeulike-article-id:11583827
* Didion, J. P., Martin, M., & Collins, F. S. (2017). Atropos: specific, sensitive, and speedy trimming of sequencing reads. <https://doi.org/10.7287/peerj.preprints.2452v4>
* Dale, R. K., Pedersen, B. S., & Quinlan, A. R. (2011). Pybedtools: A flexible Python library for manipulating genomic datasets and annotations. Bioinformatics, 27(24), 3423--3424. <https://doi.org/10.1093/bioinformatics/btr539>
* Quinlan, A. R., & Hall, I. M. (2010). BEDTools: A flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841--842. <https://doi.org/10.1093/bioinformatics/btq033>
* Tarasov, A., Vilella, A. J., Cuppen, E., Nijman, I. J., & Prins, P. (2015). Sambamba: Fast processing of NGS alignment formats. Bioinformatics, 31(12), 2032--2034. <https://doi.org/10.1093/bioinformatics/btv098>
* Heger, A. (2009). Pysam. github.com. Retrieved from <https://github.com/pysam-developers/pysam>
* Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics, 27(21), 2987--2993. <https://doi.org/10.1093/bioinformatics/btr509>
* Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., ... Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078--2079. <https://doi.org/10.1093/bioinformatics/btp352>
* Pantano, L., Estivill, X., & Martí, E. (2010). SeqBuster, a bioinformatic tool for the processing and analysis of small RNAs datasets, reveals ubiquitous miRNA modifications in human embryonic cells. Nucleic Acids Research, 38(5), e34. Retrieved from <https://www.ncbi.nlm.nih.gov/pubmed/20008100>
* Pantano, L., Friedlander, M. R., Escaramis, G., Lizano, E., Pallares-Albanell, J., Ferrer, I., ... Marti, E. (2015). Specific small-RNA signatures in the amygdala at premotor and motor stages of Parkinson's disease revealed by deep sequencing analysis. Bioinformatics (Oxford, England). <https://doi.org/10.1093/bioinformatics/btv632>

For the alignment, add what you have used:
* Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., ... Gingeras, T. R. (2013). STAR: Ultrafast universal RNA-seq aligner. Bioinformatics, 29(1), 15--21. <https://doi.org/10.1093/bioinformatics/bts635>
* Langmead, B., Trapnell, C., Pop, M., & Salzberg, S. L. (2009). Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biology, 10, R25. <https://doi.org/10.1186/gb-2009-10-3-r25>
* Kim, D., Langmead, B. & Salzberg, SL. (2016). HISAT: a fast spliced aligner with low memory requirements. Nature Methods, 12(4): 357--360. doi: 10.1038/nmeth.3317

If you used TopHat2 for alignment:
* Kim, D., Pertea, G., Trapnell, C., Pimentel, H., Kelley, R. & Salzberg SL. (2013). TopHat2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions. Genome Biology, 14(4): R36. <https://doi.org/10.1186/gb-2013-14-4-r36>
* Brueffer, C. & Saal, LH. (2016). TopHat-Recondition: A post-processor for TopHat unmapped reads. BMC Bioinformatics, 17(1):199. <https://doi.org/10.1186/s12859-016-1058-x>

If you have in the output novel miRNA discovering, add:
* Friedlander, M. R., MacKowiak, S. D., Li, N., Chen, W., & Rajewsky, N. (2012). MiRDeep2 accurately identifies known and hundreds of novel microRNA genes in seven animal clades. Nucleic Acids Research, 40(1), 37--52. <https://doi.org/10.1093/nar/gkr688>

If you have tRNA mapping output, add:
* Selitsky, S. R., & Sethupathy, P. (2015). tDRmapper: challenges and solutions to mapping, naming, and quantifying tRNA-derived RNAs from human small RNA-sequencing data. BMC Bioinformatics, 16(1), 354. <https://doi.org/10.1186/s12859-015-0800-0>

If you have miRge activated:
* Yin Lu, Alexander S. Baras, Marc K Halushka. miRge2.0: An updated tool to comprehensively analyze microRNA sequencing data. bioRxiv.org.

If you have MINTmap activated:
* Loher, P, Telonis, AG, Rigoutsos, I. MINTmap: fast and exhaustive profiling of nuclear and mitochondrial tRNA fragments from short RNA-seq data. Sci Rep. 2017;7 :41184. doi: 10.1038/srep41184. PubMed <https://www.ncbi.nlm.nih.gov/pubmed/28220888> PubMed Central PMC5318995.

### Data

* Griffiths-Jones, S. (2004). The microRNA Registry. Nucleic Acids Research, 32(Database issue), D109--11. <https://doi.org/10.1093/nar/gkh023>
* Griffiths-Jones, S. (2006). miRBase: the microRNA sequence database. Methods in Molecular Biology (Clifton, N.J.), 342, 129--38. <https://doi.org/10.1385/1-59745-123-1:129>
* Griffiths-Jones, S., Saini, H. K., Van Dongen, S., & Enright, A. J. (2008). miRBase: Tools for microRNA genomics. Nucleic Acids Research, 36(SUPPL. 1). <https://doi.org/10.1093/nar/gkm952>
* Kozomara, A., & Griffiths-Jones, S. (2011). MiRBase: Integrating microRNA annotation and deep-sequencing data. Nucleic Acids Research, 39(SUPPL. 1). <https://doi.org/10.1093/nar/gkq1027>
* Kozomara, A., & Griffiths-Jones, S. (2014). MiRBase: Annotating high confidence microRNAs using deep sequencing data. Nucleic Acids Research, 42(D1). <https://doi.org/10.1093/nar/gkt1181>
