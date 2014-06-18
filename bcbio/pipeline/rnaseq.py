import os
import bcbio.bam as bam
from bcbio.rnaseq import featureCounts, cufflinks, oncofuse, count
from bcbio.utils import get_in, safe_makedir

def detect_fusion(samples, run_parallel):
    samples = run_parallel("run_oncofuse", samples)
    return samples

def estimate_expression(samples, run_parallel):
    samples = run_parallel("generate_transcript_counts", samples)
    combined = count.combine_count_files([x[0]["count_file"] for x in samples
                                          if "count_file" in x[0]])
    gtf_file = get_in(samples[0][0], ('genome_resources', 'rnaseq',
                                      'transcripts'), None)
    annotated = count.annotate_combined_count_file(combined, gtf_file)
    samples = run_parallel("run_cufflinks", samples)
    fpkm_combined_file = os.path.splitext(combined)[0] + ".fpkm"
    to_combine = [x[0]["fpkm"] for x in samples if "fpkm" in x[0]]
    fpkm_combined = count.combine_count_files(to_combine, fpkm_combined_file)
    #fpkm_combined = cufflinks.combine_fpkm([x[0].get("fpkm_file" for x in samples]))
    for x in samples:
        x[0]["combined_counts"] = combined
        if annotated:
            x[0]["annotated_combined_counts"] = annotated
        if fpkm_combined:
            x[0]["combined_fpkm"] = fpkm_combined
    return samples

def generate_transcript_counts(data):
    """Generate counts per transcript from an alignment"""
    data["count_file"] = featureCounts.count(data)
    if get_in(data, ("config", "algorithm", "fusion_mode"), False):
        oncofuse_file = oncofuse.run(data)
        if oncofuse_file:
            data["oncofuse_file"] = oncofuse.run(data)
    return [[data]]

def run_cufflinks(data):
    """Quantitate transcript expression with Cufflinks"""
    work_bam = data["work_bam"]
    ref_file = data["sam_ref"]
    out_dir, fpkm_file = cufflinks.run(work_bam, ref_file, data)
    data["cufflinks_dir"] = out_dir
    data["fpkm"] = fpkm_file
    return [[data]]

def cufflinks_assemble(*samples):
    rnaseq_resources = samples[0][0]["genome_resources"]["rnaseq"]
    config = samples[0][0]["config"]
    dirs = samples[0][0]["dirs"]
    gtf_file = rnaseq_resources.get("transcripts", None)
    ref_file = samples[0][0]["sam_ref"]
    bam_files = [data[0]['work_bam'] for data in samples]
    num_cores = config["algorithm"].get("num_cores", 1)
    out_dir = os.path.join(dirs["work"], "assembly")
    safe_makedir(out_dir)
    merged_file = os.path.join(out_dir, "merged.bam")
    merged_file = bam.merge(bam_files, merged_file, config)
    assembly_dir = cufflinks.assemble(merged_file, ref_file, gtf_file,
                                      num_cores, out_dir)
    transcripts = [os.path.join(assembly_dir, "assembly", "transcripts.gtf")]
    merged_gtf = cufflinks.merge(transcripts, ref_file, gtf_file, num_cores)
    for data in samples:
        data[0]['assembly'] = assembly_dir
    return samples

def assemble_transcripts(run_parallel, samples):
    # assemble with all reads
    config = samples[0][0]["config"]
    if config["algorithm"].get("assemble_transcripts", False):
        samples = run_parallel("cufflinks_assemble", [samples])
    return samples
