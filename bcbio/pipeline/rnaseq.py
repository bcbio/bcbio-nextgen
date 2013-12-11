from bcbio.rnaseq import count, cufflinks

def estimate_expression(samples, run_parallel):
    samples = run_parallel("generate_transcript_counts", samples)
    samples = run_parallel("run_cufflinks", samples)
    return samples

def generate_transcript_counts(data):
    """Generate counts per transcript from an alignment"""
    data["count_file"] = count.htseq_count(data)
    return [[data]]

def run_cufflinks(data):
    """Quantitate transcript expression with Cufflinks"""
    work_bam = data["work_bam"]
    ref_file = data["sam_ref"]
    data["cufflinks_dir"] = cufflinks.run(work_bam, ref_file, data)
    return [[data]]
