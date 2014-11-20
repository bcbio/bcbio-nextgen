import os
from bcbio.rnaseq import featureCounts, cufflinks, oncofuse, count, dexseq, express
import bcbio.pipeline.datadict as dd
from bcbio.utils import filter_missing


def estimate_expression(samples, run_parallel):
    samples = run_parallel("generate_transcript_counts", samples)
    count_files = filter_missing([dd.get_count_file(x[0]) for x in samples])
    combined = count.combine_count_files(count_files)
    gtf_file = dd.get_gtf_file(samples[0][0], None)
    annotated = count.annotate_combined_count_file(combined, gtf_file)
    samples = run_parallel("run_express", samples)
    express_counts_combined = combine_express(samples, combined)

    samples = run_parallel("run_cufflinks", samples)
    #gene
    fpkm_combined_file = os.path.splitext(combined)[0] + ".fpkm"
    fpkm_files = filter_missing([dd.get_fpkm(x[0]) for x in samples])
    fpkm_combined = count.combine_count_files(fpkm_files, fpkm_combined_file)
    #isoform
    fpkm_isoform_combined_file = os.path.splitext(combined)[0] + ".isoform.fpkm"
    isoform_files = filter_missing([dd.get_fpkm_isoform(x[0]) for x in samples])
    fpkm_isoform_combined = count.combine_count_files(isoform_files,
                                                      fpkm_isoform_combined_file,
                                                      ".isoform.fpkm")
    dexseq_combined_file = os.path.splitext(combined)[0] + ".dexseq"
    to_combine_dexseq = filter_missing([dd.get_dexseq_counts(data[0]) for data in samples])
    if to_combine_dexseq:
        dexseq_combined = count.combine_count_files(to_combine_dexseq,
                                                    dexseq_combined_file, ".dexseq")
    else:
        dexseq_combined = None

    for data in dd.sample_data_iterator(samples):
        dd.set_combined_counts(data, combined)
        if annotated:
            dd.set_annotated_combined_counts(data, annotated)
        if fpkm_combined:
            dd.set_combined_fpkm(x[0], fpkm_combined)
        if fpkm_isoform_combined:
            dd.set_combined_fpkm_isoform(x[0], fpkm_combined)
        if express_counts_combined:
            dd.set_combined_eff_counts(x[0], express_counts_combined['counts'])
            dd.set_combined_tpm_counts(x[0], express_counts_combined['tpm'])
            dd.set_combined_fpkm_counts(x[0], express_counts_combined['fpkm'])
        if dexseq_combined:
            dd.set_dexseq_counts(x[0], dexseq_combined_file)
    return samples

def generate_transcript_counts(data):
    """Generate counts per transcript and per exon from an alignment"""
    data["count_file"] = featureCounts.count(data)
    if dd.get_fusion_mode(data, False):
        oncofuse_file = oncofuse.run(data)
        if oncofuse_file:
            dd.set_oncofuse_file(data, oncofuse_file)
    if dd.get_dexseq_gff(data, None):
        data = dd.set_dexseq_counts(data, dexseq.bcbio_run(data))
    # if RSEM was run, stick the transcriptome BAM file into the datadict
    if dd.get_aligner(data).lower() == "star" and dd.get_rsem(data):
        base, ext = os.path.splitext(dd.get_work_bam(data))
        data = dd.set_transcriptome_bam(data, base + ".transcriptome" + ext)
    return [[data]]

def run_express(data):
    """Quantitative isoform expression by  express"""
    out_files = express.run(data)
    if out_files:
        data['eff_counts'], data['tpm_counts'], data['fpkm_counts'] = out_files
    return [[data]]

def combine_express(samples, combined):
    """Combine tpm, effective counts and fpkm from express results"""
    to_combine = [x[0]["eff_counts"] for x in samples if "eff_counts" in x[0]]
    if len(to_combine) > 0:
        eff_counts_combined_file = os.path.splitext(combined)[0] + "_eff.counts"
        eff_counts_combined = count.combine_count_files(to_combine, eff_counts_combined_file)
        to_combine = [x[0]["tpm_counts"] for x in samples if "tpm_counts" in x[0]]
        tpm_counts_combined_file = os.path.splitext(combined)[0] + ".tpm"
        tpm_counts_combined = count.combine_count_files(to_combine, tpm_counts_combined_file)
        to_combine = [x[0]["fpkm_counts"] for x in samples if "fpkm_counts" in x[0]]
        fpkm_counts_combined_file = os.path.splitext(combined)[0] + ".fpkm"
        fpkm_counts_combined = count.combine_count_files(to_combine, fpkm_counts_combined_file)
        return {'counts': eff_counts_combined, 'tpm': tpm_counts_combined,
                'fpkm': fpkm_counts_combined}
    return None

def run_cufflinks(data):
    """Quantitate transcript expression with Cufflinks"""
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
    data = dd.set_assembled_gtf(data, assembled_gtf)
    return [[data]]

def cufflinks_merge(*samples):
    to_merge = filter_missing([dd.get_assembled_gtf(data) for data in
                            dd.sample_data_iterator(samples)])
    data = samples[0][0]
    bam_file = dd.get_work_bam(data)
    ref_file = dd.get_sam_ref(data)
    gtf_file = dd.get_gtf_file(data)
    out_dir = os.path.join(dd.get_work_dir(data), "assembly")
    num_cores = dd.get_num_cores(data)
    merged_gtf = cufflinks.merge(to_merge, ref_file, gtf_file, num_cores, samples[0][0])
    for data in dd.sample_data_iterator(samples):
        dd.set_assembled_gtf(data, merged_gtf)
    return samples

def assemble_transcripts(run_parallel, samples):
    """
    assembly strategy rationale implemented as suggested in
    http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html

    run Cufflinks in without a reference GTF for each individual sample
    merge the assemblies with Cuffmerge using a reference GTF
    """
    if dd.get_assemble_transcripts(samples[0][0]):
        samples = run_parallel("cufflinks_assemble", samples)
        samples = run_parallel("cufflinks_merge", [samples])
    return samples
