import os
import sys
from bcbio.rnaseq import (featureCounts, cufflinks, oncofuse, count, dexseq,
                          express, variation, stringtie, sailfish, spikein, pizzly, ericscript,
                          kallisto, salmon, singlecellexperiment, arriba)
from bcbio.ngsalign import bowtie2, alignprep
from bcbio.variation import effects, joint, multi, population, vardict
import bcbio.pipeline.datadict as dd
from bcbio.utils import filter_missing, flatten, to_single_data, file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger

def fast_rnaseq(samples, run_parallel):
    samples = run_parallel("run_salmon_index", [samples])
    samples = run_parallel("run_salmon_reads", samples)
    samples = run_parallel("run_counts_spikein", samples)
    samples = spikein.combine_spikein(samples)
    return samples

def singlecell_rnaseq(samples, run_parallel):
    quantifier = dd.get_in_samples(samples, dd.get_singlecell_quantifier)
    quantifier = quantifier.lower()
    samples = run_parallel("run_umi_transform", samples)
    demultiplexed = run_parallel("demultiplex_samples", samples)
    # break demultiplixed lanes into their own samples
    samples = []
    for lane in demultiplexed:
        for index in lane:
            samples.append([index])
    samples = run_parallel("run_filter_barcodes", samples)
    samples = run_parallel("run_barcode_histogram", samples)
    if quantifier == "rapmap":
        samples = run_parallel("run_rapmap_index", [samples])
        samples = run_parallel("run_rapmap_align", samples)
        samples = run_parallel("run_tagcount", samples)
        samples = run_parallel("run_concatenate_sparse_counts", [samples])
    elif quantifier == "kallisto":
        samples = run_parallel("run_kallisto_singlecell", samples)
    else:
        logger.error(("%s is not supported for singlecell RNA-seq "
                      "quantification." % quantifier))
        sys.exit(1)
    samples = scrnaseq_concatenate_metadata(samples)
    singlecellexperiment.make_scrnaseq_object(samples)
    return samples

def scrnaseq_concatenate_metadata(samples):
    """
    Create file same dimension than mtx.colnames
    with metadata and sample name to help in the
    creation of the SC object.
    """
    barcodes = {}
    counts =  ""
    metadata = {}
    has_sample_barcodes = False
    for sample in dd.sample_data_iterator(samples):
        if dd.get_sample_barcodes(sample):
            has_sample_barcodes = True
            with open(dd.get_sample_barcodes(sample)) as inh:
                for line in inh:
                    cols = line.strip().split(",")
                    if len(cols) == 1:
                        # Assign sample name in case of missing in barcodes
                        cols.append("NaN")
                    barcodes[(dd.get_sample_name(sample), cols[0])] = cols[1:]
        else:
            barcodes[(dd.get_sample_name(sample), "NaN")] = [dd.get_sample_name(sample), "NaN"]

        counts = dd.get_combined_counts(sample)
        meta = map(str, list(sample["metadata"].values()))
        meta_cols = list(sample["metadata"].keys())
        meta = ["NaN" if not v else v for v in meta]
        metadata[dd.get_sample_name(sample)] = meta

    metadata_fn = counts + ".metadata"
    if file_exists(metadata_fn):
        return samples
    with file_transaction(metadata_fn) as tx_metadata_fn:
        with open(tx_metadata_fn, 'w') as outh:
            outh.write(",".join(["sample"] + meta_cols) + '\n')
            with open(counts + ".colnames") as inh:
                for line in inh:
                    sample = line.split(":")[0]
                    if has_sample_barcodes:
                        barcode = sample.split("-")[1]
                    else:
                        barcode = "NaN"
                    outh.write(",".join(barcodes[(sample, barcode)] + metadata[sample]) + '\n')
    return samples

def rnaseq_variant_calling(samples, run_parallel):
    """
    run RNA-seq variant calling using GATK
    """
    samples = run_parallel("run_rnaseq_variant_calling", samples)
    variantcaller = dd.get_variantcaller(to_single_data(samples[0]))
    if variantcaller and ("gatk-haplotype" in variantcaller):
        out = []
        for d in joint.square_off(samples, run_parallel):
            out.extend([[to_single_data(xs)] for xs in multi.split_variants_by_sample(to_single_data(d))])
        samples = out
    if variantcaller:
        samples = run_parallel("run_rnaseq_ann_filter", samples)
    if variantcaller and ("gatk-haplotype" in variantcaller):
        out = []
        for data in (to_single_data(xs) for xs in samples):
            if "variants" not in data:
                data["variants"] = []
            data["variants"].append({"variantcaller": "gatk-haplotype", "vcf": data["vrn_file_orig"],
                                     "population": {"vcf": data["vrn_file"]}})
            data["vrn_file"] = data.pop("vrn_file_orig")
            out.append([data])
        samples = out
    return samples

def run_rnaseq_variant_calling(data):
    """
    run RNA-seq variant calling, variation file is stored in `vrn_file`
    in the datadict
    """
    variantcaller = dd.get_variantcaller(data)
    if isinstance(variantcaller, list) and len(variantcaller) > 1:
        logger.error("Only one variantcaller can be run for RNA-seq at "
                     "this time. Post an issue here "
                     "(https://github.com/bcbio/bcbio-nextgen/issues) "
                     "if this is something you need to do.")
        sys.exit(1)

    if variantcaller:
        if "gatk-haplotype" in variantcaller:
            data = variation.rnaseq_gatk_variant_calling(data)
        if vardict.get_vardict_command(data):
            data = variation.rnaseq_vardict_variant_calling(data)
        vrn_file = dd.get_vrn_file(data)
    return [[data]]

def run_rnaseq_ann_filter(data):
    """Run RNA-seq annotation and filtering.
    """
    data = to_single_data(data)
    if dd.get_vrn_file(data):
        eff_file = effects.add_to_vcf(dd.get_vrn_file(data), data)[0]
        if eff_file:
            data = dd.set_vrn_file(data, eff_file)
        ann_file = population.run_vcfanno(dd.get_vrn_file(data), data)
        if ann_file:
            data = dd.set_vrn_file(data, ann_file)
    variantcaller = dd.get_variantcaller(data)
    if variantcaller and ("gatk-haplotype" in variantcaller):
        filter_file = variation.gatk_filter_rnaseq(dd.get_vrn_file(data), data)
        data = dd.set_vrn_file(data, filter_file)
    # remove variants close to splice junctions
    vrn_file = dd.get_vrn_file(data)
    vrn_file = variation.filter_junction_variants(vrn_file, data)
    data = dd.set_vrn_file(data, vrn_file)
    return [[data]]

def quantitate(data):
    """CWL target for quantitation.

    XXX Needs to be split and parallelized by expression caller, with merging
    of multiple calls.
    """
    data = to_single_data(to_single_data(data))
    data = generate_transcript_counts(data)[0][0]
    data["quant"] = {}
    if "sailfish" in dd.get_expression_caller(data):
        data = to_single_data(sailfish.run_sailfish(data)[0])
        data["quant"]["tsv"] = data["sailfish"]
        data["quant"]["hdf5"] = os.path.join(os.path.dirname(data["sailfish"]), "abundance.h5")
    if ("kallisto" in dd.get_expression_caller(data) or "pizzly" in dd.get_fusion_caller(data, [])):
        data = to_single_data(kallisto.run_kallisto_rnaseq(data)[0])
        data["quant"]["tsv"] = os.path.join(data["kallisto_quant"], "abundance.tsv")
        data["quant"]["hdf5"] = os.path.join(data["kallisto_quant"], "abundance.h5")
    if (os.path.exists(os.path.join(data["kallisto_quant"], "fusion.txt"))):
        data["quant"]["fusion"] = os.path.join(data["kallisto_quant"], "fusion.txt")
    else:
        data["quant"]["fusion"] = None
    if "salmon" in dd.get_expression_caller(data):
        if dd.get_quantify_genome_alignments(data): 
            if dd.get_aligner(data).lower() != "star":
                if dd.get_genome_build(data) == "hg38":
                    logger.warning("Whole genome alignment-based Salmon quantification is "
                         "only supported for the STAR aligner. Since this is hg38 we will fall "
                         "back to the decoy method")
                    data = to_single_data(salmon.run_salmon_decoy(data)[0])
                else:
                    logger.warning(
                         "Whole genome alignment-based Salmon quantification is "
                         "only supported for the STAR aligner. Falling back to the "
                         "transcriptome-only method.")
                    data = to_single_data(salmon.run_salmon_reads(data)[0])
            else:
                data = to_single_data(salmon.run_salmon_bam(data)[0])
        else:
            data = to_single_data(salmon.run_salmon_reads(data)[0])
        data["quant"]["tsv"] = data["salmon"]
        data["quant"]["hdf5"] = os.path.join(os.path.dirname(data["salmon"]), "abundance.h5")
    return [[data]]

def quantitate_expression_parallel(samples, run_parallel):
    """
    quantitate expression, all programs run here should be multithreaded to
    take advantage of the threaded run_parallel environment
    """
    data = samples[0][0]
    samples = run_parallel("generate_transcript_counts", samples)
    if "cufflinks" in dd.get_expression_caller(data):
        samples = run_parallel("run_cufflinks", samples)
    if "stringtie" in dd.get_expression_caller(data):
        samples = run_parallel("run_stringtie_expression", samples)
    if ("kallisto" in dd.get_expression_caller(data) or
        dd.get_fusion_mode(data) or
        "pizzly" in dd.get_fusion_caller(data, [])):
        samples = run_parallel("run_kallisto_index", [samples])
        samples = run_parallel("run_kallisto_rnaseq", samples)
    if "sailfish" in dd.get_expression_caller(data):
        samples = run_parallel("run_sailfish_index", [samples])
        samples = run_parallel("run_sailfish", samples)

    # always run salmon
    if dd.get_quantify_genome_alignments(data):
        if dd.get_aligner(data).lower() != "star":
            if dd.get_genome_build(data) == "hg38":
                logger.warning("Whole genome alignment-based Salmon quantification is "
                   "only supported for the STAR aligner. Since this is hg38 we will fall "
                   "back to the decoy method")
                samples = run_parallel("run_salmon_decoy", samples)
            else:
                logger.warning(
                   "Whole genome alignment-based Salmon quantification is "
                   "only supported for the STAR aligner. Falling back to the "
                   "transcriptome-only method.")
                samples = run_parallel("run_salmon_reads", samples)
        else:
            samples = run_parallel("run_salmon_bam", samples)
    else:
        samples = run_parallel("run_salmon_reads", samples)

    samples = run_parallel("detect_fusions", samples)
    return samples

def detect_fusions(data):
    data = to_single_data(data)
    # support the old style of fusion mode calling
    if dd.get_fusion_mode(data, False):
        data = dd.set_fusion_caller(data, ["oncofuse", "pizzly"])
        logger.warning("``fusion_mode`` is deprecated in favor of turning on "
                       "callers with ``fusion_caller``. It will run pizzly and "
                       "oncofuse for now, but will eventually have support "
                       "dropped.")
    fusion_caller = dd.get_fusion_caller(data, [])
    if "oncofuse" in fusion_caller:
        oncofuse_file = oncofuse.run(data)
        if oncofuse_file:
            data = dd.set_oncofuse_file(data, oncofuse_file)
    if "pizzly" in fusion_caller:
        pizzly_dir = pizzly.run_pizzly(data)
        if pizzly_dir:
            data = dd.set_pizzly_dir(data, pizzly_dir)
            data["fusion"] = {"fasta": os.path.join(pizzly_dir, "%s.fusions.fasta" % dd.get_sample_name(data)),
                              "json": os.path.join(pizzly_dir, "%s.json" % dd.get_sample_name(data))}
    if "ericscript" in fusion_caller:
        ericscript_dir = ericscript.run(data)
    if "arriba" in fusion_caller:
        data = arriba.run_arriba(data)
    return [[data]]

def quantitate_expression_noparallel(samples, run_parallel):
    """
    run transcript quantitation for algorithms that don't run in parallel
    """
    data = samples[0][0]
    if "express" in dd.get_expression_caller(data):
        samples = run_parallel("run_express", samples)
    if "dexseq" in dd.get_expression_caller(data):
        samples = run_parallel("run_dexseq", samples)
    return samples

def generate_transcript_counts(data):
    """Generate counts per transcript and per exon from an alignment"""
    data["count_file"] = featureCounts.count(data)

    if dd.get_fusion_mode(data, False) and not dd.get_fusion_caller(data):
        oncofuse_file = oncofuse.run(data)
        if oncofuse_file:
            data = dd.set_oncofuse_file(data, oncofuse_file)

    if dd.get_transcriptome_align(data):
        # to create a disambiguated transcriptome file realign with bowtie2
        if dd.get_disambiguate(data):
            logger.info("Aligning to the transcriptome with bowtie2 using the "
                        "disambiguated reads.")
            bam_path = data["work_bam"]
            fastq_paths = alignprep._bgzip_from_bam(bam_path, data["dirs"], data, is_retry=False, output_infix='-transcriptome')
            if len(fastq_paths) == 2:
                file1, file2 = fastq_paths
            else:
                file1, file2 = fastq_paths[0], None
            ref_file = dd.get_ref_file(data)
            data = bowtie2.align_transcriptome(file1, file2, ref_file, data)
        else:
            file1, file2 = dd.get_input_sequence_files(data)
        if not dd.get_transcriptome_bam(data):
            ref_file = dd.get_ref_file(data)
            logger.info("Transcriptome alignment was flagged to run, but the "
                        "transcriptome BAM file was not found. Aligning to the "
                        "transcriptome with bowtie2.")
            data = bowtie2.align_transcriptome(file1, file2, ref_file, data)
    data = spikein.counts_spikein(data)
    return [[data]]

def run_stringtie_expression(data):
    """Calculate transcript and gene level FPKM with Stringtie"""
    data = stringtie.run_stringtie_expression(data)
    return [[data]]

def run_dexseq(data):
    """Quantitate exon-level counts with DEXSeq"""
    if dd.get_dexseq_gff(data, None):
        data = dexseq.bcbio_run(data)
    return [[data]]

def run_express(data):
    """Quantitative isoform expression by eXpress"""
    data = express.run(data)
    return [[data]]

def combine_express(samples, combined):
    """Combine tpm, effective counts and fpkm from express results"""
    if not combined:
        return None
    to_combine = [dd.get_express_counts(x) for x in
                  dd.sample_data_iterator(samples) if dd.get_express_counts(x)]
    gtf_file = dd.get_gtf_file(samples[0][0])
    isoform_to_gene_file = os.path.join(os.path.dirname(combined), "isoform_to_gene.txt")
    isoform_to_gene_file = express.isoform_to_gene_name(
        gtf_file, isoform_to_gene_file, next(dd.sample_data_iterator(samples)))
    if len(to_combine) > 0:
        eff_counts_combined_file = os.path.splitext(combined)[0] + ".isoform.express_counts"
        eff_counts_combined = count.combine_count_files(to_combine, eff_counts_combined_file, ext=".counts")
        to_combine = [dd.get_express_tpm(x) for x in
                      dd.sample_data_iterator(samples) if dd.get_express_tpm(x)]
        tpm_counts_combined_file = os.path.splitext(combined)[0] + ".isoform.express_tpm"
        tpm_counts_combined = count.combine_count_files(to_combine, tpm_counts_combined_file)
        to_combine = [dd.get_express_fpkm(x) for x in dd.sample_data_iterator(samples)
                      if dd.get_express_fpkm(x)]
        fpkm_counts_combined_file = os.path.splitext(combined)[0] + ".isoform.express_fpkm"
        fpkm_counts_combined = count.combine_count_files(to_combine, fpkm_counts_combined_file, ext=".fpkm")
        return {'counts': eff_counts_combined, 'tpm': tpm_counts_combined,
                'fpkm': fpkm_counts_combined, 'isoform_to_gene': isoform_to_gene_file}
    return {}

def run_cufflinks(data):
    """Quantitate transcript expression with Cufflinks"""
    if "cufflinks" in dd.get_tools_off(data):
        return [[data]]
    work_bam = dd.get_work_bam(data)
    ref_file = dd.get_sam_ref(data)
    out_dir, fpkm_file, fpkm_isoform_file = cufflinks.run(work_bam, ref_file, data)
    data = dd.set_cufflinks_dir(data, out_dir)
    data = dd.set_fpkm(data, fpkm_file)
    data = dd.set_fpkm_isoform(data, fpkm_isoform_file)
    return [[data]]

def cufflinks_assemble(data):
    bam_file = dd.get_work_bam(data)
    ref_file = dd.get_sam_ref(data)
    out_dir = os.path.join(dd.get_work_dir(data), "assembly")
    num_cores = dd.get_num_cores(data)
    assembled_gtf = cufflinks.assemble(bam_file, ref_file, num_cores, out_dir, data)
    dd.get_assembled_gtf(data).append(assembled_gtf)
    return [[data]]

def cufflinks_merge(*samples):
    to_merge = set(filter_missing(flatten([dd.get_assembled_gtf(data) for data in
                                           dd.sample_data_iterator(samples)])))
    data = samples[0][0]
    ref_file = dd.get_sam_ref(data)
    gtf_file = dd.get_gtf_file(data)
    num_cores = dd.get_num_cores(data)
    merged_gtf = cufflinks.merge(to_merge, ref_file, gtf_file, num_cores,
                                 samples[0][0])
    updated_samples = []
    for data in dd.sample_data_iterator(samples):
        data = dd.set_merged_gtf(data, merged_gtf)
        updated_samples.append([data])
    return updated_samples

def stringtie_merge(*samples):
    to_merge = set(filter_missing(flatten([dd.get_assembled_gtf(data) for data in
                                       dd.sample_data_iterator(samples)])))
    data = samples[0][0]
    ref_file = dd.get_sam_ref(data)
    gtf_file = dd.get_gtf_file(data)
    num_cores = dd.get_num_cores(data)
    merged_gtf = stringtie.merge(to_merge, ref_file, gtf_file, num_cores, data)
    updated_samples = []
    for data in dd.sample_data_iterator(samples):
        data = dd.set_merged_gtf(data, merged_gtf)
        updated_samples.append([data])
    return updated_samples

def assemble_transcripts(run_parallel, samples):
    """
    assembly strategy rationale implemented as suggested in
    http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html

    run Cufflinks in without a reference GTF for each individual sample
    merge the assemblies with Cuffmerge using a reference GTF
    """
    assembler = dd.get_in_samples(samples, dd.get_transcript_assembler)
    data = samples[0][0]
    if assembler:
        if "cufflinks" in assembler:
            samples = run_parallel("cufflinks_assemble", samples)
        if "stringtie" in assembler:
            samples = run_parallel("run_stringtie_expression", samples)
        if "stringtie" in assembler and stringtie.supports_merge(data):
            samples = run_parallel("stringtie_merge", [samples])
        else:
            samples = run_parallel("cufflinks_merge", [samples])
    return samples

def combine_files(samples):
    """
    after quantitation, combine the counts/FPKM/TPM/etc into a single table with
    all samples
    """
    data = samples[0][0]
    # prefer the supplied transcriptome gtf file
    gtf_file = dd.get_transcriptome_gtf(data, None)
    if not gtf_file:
        gtf_file = dd.get_gtf_file(data, None)
    dexseq_gff = dd.get_dexseq_gff(data)

    # combine featureCount files
    count_files = filter_missing([dd.get_count_file(x[0]) for x in samples])
    combined = count.combine_count_files(count_files, ext=".counts")
    annotated = count.annotate_combined_count_file(combined, gtf_file)

    # add tx2gene file
    tx2gene_file = os.path.join(dd.get_work_dir(data), "annotation", "tx2gene.csv")
    if gtf_file:
        tx2gene_file = sailfish.create_combined_tx2gene(data)

    # combine eXpress files
    express_counts_combined = combine_express(samples, combined)

    # combine Cufflinks files
    fpkm_files = filter_missing([dd.get_fpkm(x[0]) for x in samples])
    if fpkm_files and combined:
        fpkm_combined_file = os.path.splitext(combined)[0] + ".fpkm"
        fpkm_combined = count.combine_count_files(fpkm_files, fpkm_combined_file)
    else:
        fpkm_combined = None
    isoform_files = filter_missing([dd.get_fpkm_isoform(x[0]) for x in samples])
    if isoform_files and combined:
        fpkm_isoform_combined_file = os.path.splitext(combined)[0] + ".isoform.fpkm"
        fpkm_isoform_combined = count.combine_count_files(isoform_files,
                                                          fpkm_isoform_combined_file,
                                                          ".isoform.fpkm")
    else:
        fpkm_isoform_combined = None
    # combine DEXseq files
    to_combine_dexseq = filter_missing([dd.get_dexseq_counts(data[0]) for data
                                        in samples])
    if to_combine_dexseq and combined:
        dexseq_combined_file = os.path.splitext(combined)[0] + ".dexseq"
        dexseq_combined = count.combine_count_files(to_combine_dexseq,
                                                    dexseq_combined_file, ".dexseq")
        if dexseq_combined:
            dexseq.create_dexseq_annotation(dexseq_gff, dexseq_combined)
    else:
        dexseq_combined = None
    samples = spikein.combine_spikein(samples)
    updated_samples = []
    for data in dd.sample_data_iterator(samples):
        if combined:
            data = dd.set_combined_counts(data, combined)
        if annotated:
            data = dd.set_annotated_combined_counts(data, annotated)
        if fpkm_combined:
            data = dd.set_combined_fpkm(data, fpkm_combined)
        if fpkm_isoform_combined:
            data = dd.set_combined_fpkm_isoform(data, fpkm_isoform_combined)
        if express_counts_combined:
            data = dd.set_express_counts(data, express_counts_combined['counts'])
            data = dd.set_express_tpm(data, express_counts_combined['tpm'])
            data = dd.set_express_fpkm(data, express_counts_combined['fpkm'])
            data = dd.set_isoform_to_gene(data, express_counts_combined['isoform_to_gene'])
        if dexseq_combined:
            data = dd.set_dexseq_counts(data, dexseq_combined_file)
        if gtf_file:
            data = dd.set_tx2gene(data, tx2gene_file)
        updated_samples.append([data])
    return updated_samples
