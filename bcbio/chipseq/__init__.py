from bcbio.ngsalign.bowtie2 import filter_multimappers
import bcbio.pipeline.datadict as dd
from bcbio.log import logger

def clean_chipseq_alignment(data):
    aligner = dd.get_aligner(data)
    data["raw_bam"] = dd.get_work_bam(data)
    if aligner:
        assert aligner == "bowtie2", "ChIP-seq only supported for bowtie2."
        unique_bam = filter_multimappers(dd.get_work_bam(data), data)
        data["work_bam"] = unique_bam
    else:
        logger.info("Warning: When BAM file is given as input, bcbio skips multimappers removal."
                    "If BAM is not cleaned for peak calling, can result in downstream errors.")
    return [[data]]
