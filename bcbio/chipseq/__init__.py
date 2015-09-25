from bcbio.ngsalign.bowtie2 import filter_multimappers
import bcbio.pipeline.datadict as dd

def clean_chipseq_alignment(data):
    aligner = dd.get_aligner(data)
    assert aligner == "bowtie2", "ChIP-seq only supported for bowtie2."
    if aligner == "bowtie2":
        unique_bam = filter_multimappers(dd.get_work_bam(data), data)
        data["work_bam"] = unique_bam
        return [[data]]
    return [[data]]
