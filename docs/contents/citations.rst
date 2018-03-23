Citations
---------

https://github.com/bcbio/bcbio-nextgen

small RNA-seq
=============

Data was analyzed with bcbio-nextgen (https://github.com/bcbio/bcbio-nextgen)
using piDNA to detect the adapter, cutadapt to remove it, STAR/bowtie to align against
the genome and seqcluster to detect small RNA transcripts. miRNAs were detected using
miraligner tool with miRBase as the reference miRNA database. tRNA profiles were
detected using tdrmapper tool. mirdeep2 was used for discovery of novel miRNAs. FastQC
was used for QC metrics and multiqc for reporting.

Download BIB format: https://github.com/bcbio/bcbio-nextgen/tree/master/docs/contents/misc/bcbio-smallrna.bib

Tools
~~~~~

* Tsuji J, Weng Z. (2016) DNApi: A De Novo Adapter Prediction Algorithm for Small
  RNA Sequencing Data. 11(10):e0164228. http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0164228

* Andrews, S. (2010). FastQC: A quality control tool for high throughput sequence data. Bioinformatics. doi:citeulike-article-id:11583827

* Didion, J. P., Martin, M., & Collins, F. S. (2017). Atropos: specific, sensitive, and speedy trimming of sequencing reads. http://doi.org/10.7287/peerj.preprints.2452v4

* Dale, R. K., Pedersen, B. S., & Quinlan, A. R. (2011). Pybedtools: A flexible Python library for manipulating genomic datasets and annotations. Bioinformatics, 27(24), 3423–3424. doi:10.1093/bioinformatics/btr539

* Quinlan, A. R., & Hall, I. M. (2010). BEDTools: A flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841–842. doi:10.1093/bioinformatics/btq033

* Tarasov, A., Vilella, A. J., Cuppen, E., Nijman, I. J., & Prins, P. (2015). Sambamba: Fast processing of NGS alignment formats. Bioinformatics, 31(12), 2032–2034. doi:10.1093/bioinformatics/btv098

* Heger, A. (2009). Pysam. github.com. Retrieved from https://github.com/pysam-developers/pysam

* Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics, 27(21), 2987–2993. doi:10.1093/bioinformatics/btr509

* Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., … Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. doi:10.1093/bioinformatics/btp352

* Pantano, L., Estivill, X., & Martí, E. (2010). SeqBuster, a bioinformatic tool for the processing and analysis of small RNAs datasets, reveals ubiquitous miRNA modifications in human embryonic cells. Nucleic Acids Research, 38(5), e34. Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/20008100

* Pantano, L., Friedlander, M. R., Escaramis, G., Lizano, E., Pallares-Albanell, J., Ferrer, I., … Marti, E. (2015). Specific small-RNA signatures in the amygdala at premotor and motor stages of Parkinson’s disease revealed by deep sequencing analysis. Bioinformatics (Oxford, England). doi:10.1093/bioinformatics/btv632


For the alignment, add what you have used:

* Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., … Gingeras, T. R. (2013). STAR: Ultrafast universal RNA-seq aligner. Bioinformatics, 29(1), 15–21. doi:10.1093/bioinformatics/bts635

* Langmead, B., Trapnell, C., Pop, M., & Salzberg, S. L. (2009). Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biology, 10, R25. doi:10.1186/gb-2009-10-3-r25

* Kim, D., Langmead, B. & Salzberg, SL. (2016). HISAT: a fast spliced aligner with low memory requirements. Nature Methods, 12(4): 357–360. doi: 10.1038/nmeth.3317


If you used TopHat2 for alignment:

* Kim, D., Pertea, G., Trapnell, C., Pimentel, H., Kelley, R. & Salzberg SL. (2013). TopHat2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions. Genome Biology, 14(4): R36. doi: 10.1186/gb-2013-14-4-r36

* Brueffer, C. & Saal, LH. (2016). TopHat-Recondition: A post-processor for TopHat unmapped reads. BMC Bioinformatics, 17(1):199. doi: 10.1186/s12859-016-1058-x


If you have in the output novel miRNA discovering, add: 

* Friedlander, M. R., MacKowiak, S. D., Li, N., Chen, W., & Rajewsky, N. (2012). MiRDeep2 accurately identifies known and hundreds of novel microRNA genes in seven animal clades. Nucleic Acids Research, 40(1), 37–52. doi:10.1093/nar/gkr688

If you have tRNA mapping output, add:

* Selitsky, S. R., & Sethupathy, P. (2015). tDRmapper: challenges and solutions to mapping, naming, and quantifying tRNA-derived RNAs from human small RNA-sequencing data. BMC Bioinformatics, 16(1), 354. doi:10.1186/s12859-015-0800-0

If you have miRge activated:

* Yin Lu, Alexander S. Baras, Marc K Halushka. miRge2.0: An updated tool to comprehensively analyze microRNA sequencing data. bioRxiv.org.

If you have MINTmap activated:

* Loher, P, Telonis, AG, Rigoutsos, I. MINTmap: fast and exhaustive profiling of nuclear and mitochondrial tRNA fragments from short RNA-seq data. Sci Rep. 2017;7 :41184. doi: 10.1038/srep41184. PubMed PMID:28220888 PubMed Central PMC5318995.

Data
~~~~

* Griffiths-Jones, S. (2004). The microRNA Registry. Nucleic Acids Research, 32(Database issue), D109–11. doi:10.1093/nar/gkh023

* Griffiths-Jones, S. (2006). miRBase: the microRNA sequence database. Methods in Molecular Biology (Clifton, N.J.), 342, 129–38. doi:10.1385/1-59745-123-1:129

* Griffiths-Jones, S., Saini, H. K., Van Dongen, S., & Enright, A. J. (2008). miRBase: Tools for microRNA genomics. Nucleic Acids Research, 36(SUPPL. 1). doi:10.1093/nar/gkm952

* Kozomara, A., & Griffiths-Jones, S. (2011). MiRBase: Integrating microRNA annotation and deep-sequencing data. Nucleic Acids Research, 39(SUPPL. 1). doi:10.1093/nar/gkq1027

* Kozomara, A., & Griffiths-Jones, S. (2014). MiRBase: Annotating high confidence microRNAs using deep sequencing data. Nucleic Acids Research, 42(D1). doi:10.1093/nar/gkt1181
