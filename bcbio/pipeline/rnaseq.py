import os
from bcbio.rnaseq import featureCounts, cufflinks, oncofuse, count, dexseq
import bcbio.pipeline.datadict as dd

def estimate_expression(samples, run_parallel):
    samples = run_parallel("generate_transcript_counts", samples)
    combined = count.combine_count_files([x[0]["count_file"] for x in samples
                                          if "count_file" in x[0]])
    gtf_file = dd.get_gtf_file(samples[0][0], None)
    annotated = count.annotate_combined_count_file(combined, gtf_file)
    samples = run_parallel("run_cufflinks", samples)
    #gene
    fpkm_combined_file = os.path.splitext(combined)[0] + ".fpkm"
    to_combine = [x[0]["fpkm"] for x in samples if "fpkm" in x[0]]
    fpkm_combined = count.combine_count_files(to_combine, fpkm_combined_file)
    #isoform
    fpkm_isoform_combined_file = os.path.splitext(combined)[0] + ".isoform.fpkm"
    to_combine_isoform = [x[0]["fpkm_isoform"] for x in samples if "fpkm_isoform" in x[0]]
    fpkm_isoform_combined = count.combine_count_files(to_combine_isoform,
                                                      fpkm_isoform_combined_file,
                                                      ".isoform.fpkm")
    dexseq_combined_file = os.path.splitext(combined)[0] + ".dexseq"
    to_combine_dexseq = [dd.get_dexseq_counts(data[0]) for data in samples]
    to_combine_dexseq = filter(lambda x: x, to_combine_dexseq)
    if to_combine_dexseq:
        dexseq_combined = count.combine_count_files(to_combine_dexseq,
                                                    dexseq_combined_file, ".dexseq")
    else:
        dexseq_combined = None

    for x in samples:
        x[0]["combined_counts"] = combined
        if annotated:
            x[0]["annotated_combined_counts"] = annotated
        if fpkm_combined:
            x[0]["combined_fpkm"] = fpkm_combined
        if fpkm_isoform_combined:
            x[0]["combined_fpkm_isoform"] = fpkm_isoform_combined
        if dexseq_combined:
            x[0] = dd.set_dexseq_counts(x[0], dexseq_combined_file)

    return samples

def generate_transcript_counts(data):
    """Generate counts per transcript and per exon from an alignment"""
    data["count_file"] = featureCounts.count(data)
    if dd.get_fusion_mode(data, False):
        oncofuse_file = oncofuse.run(data)
        if oncofuse_file:
            data["oncofuse_file"] = oncofuse.run(data)
    if dd.get_dexseq_gff(data, None):
        data = dd.set_dexseq_counts(data, dexseq.bcbio_run(data))
    return [[data]]

def run_cufflinks(data):
    """Quantitate transcript expression with Cufflinks"""
    work_bam = data["work_bam"]
    ref_file = data["sam_ref"]
    out_dir, fpkm_file, fpkm_isoform_file = cufflinks.run(work_bam, ref_file, data)
    data["cufflinks_dir"] = out_dir
    data["fpkm"] = fpkm_file
    data["fpkm_isoform"] = fpkm_isoform_file
    return [[data]]

def cufflinks_assemble(data):
    config = data["config"]
    dirs = data["dirs"]
    bam_file = data["work_bam"]
    ref_file = data["sam_ref"]
    out_dir = os.path.join(dirs["work"], "assembly")
    num_cores = config["algorithm"].get("num_cores", 1)
    assembled_gtf = cufflinks.assemble(bam_file, ref_file, num_cores, out_dir, data)
    data["assembled_gtf"] = assembled_gtf
    return [[data]]

def cufflinks_merge(*samples):
    rnaseq_resources = samples[0][0]["genome_resources"]["rnaseq"]
    config = samples[0][0]["config"]
    dirs = samples[0][0]["dirs"]
    bam_file = samples[0][0]["work_bam"]
    ref_file = samples[0][0]["sam_ref"]
    gtf_file = rnaseq_resources.get("transcripts", None)
    out_dir = os.path.join(dirs["work"], "assembly")
    num_cores = config["algorithm"].get("num_cores", 1)
    to_merge = [data[0]["assembled_gtf"] for data in samples if
                "assembled_gtf" in data[0]]
    merged_gtf = cufflinks.merge(to_merge, ref_file, gtf_file, num_cores, samples[0][0])
    for data in samples:
        data[0]['assembled_gtf'] = merged_gtf
    return samples

def assemble_transcripts(run_parallel, samples):
    """
    assembly strategy rationale implemented as suggested in
    http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html

    run Cufflinks in without a reference GTF for each individual sample
    merge the assemblies with Cuffmerge using a reference GTF
    """
    config = samples[0][0]["config"]
    if config["algorithm"].get("assemble_transcripts", False):
        samples = run_parallel("cufflinks_assemble", samples)
        samples = run_parallel("cufflinks_merge", [samples])
    return samples
