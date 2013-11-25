from bcbio.ngsalign.bowtie2 import filter_multimappers


def clean_chipseq_alignment(data):
    config = data["config"]
    aligner = config["algorithm"].get("aligner", None)
    assert aligner == "bowtie2", "ChIP-seq only supported for bowtie2."
    if aligner == "bowtie2":
        work_bam = filter_multimappers(data["work_bam"])
        data["work_bam"] = work_bam
        return [[data]]
