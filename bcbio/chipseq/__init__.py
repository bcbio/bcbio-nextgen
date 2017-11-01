import sys
import bcbio.pipeline.datadict as dd
from bcbio.ngsalign import bowtie2, bwa
from bcbio.log import logger

def clean_chipseq_alignment(data):
    aligner = dd.get_aligner(data)
    data["raw_bam"] = dd.get_work_bam(data)
    if aligner:
        if aligner == "bowtie2":
            filterer = bowtie2.filter_multimappers
        elif aligner == "bwa":
            filterer = bwa.filter_multimappers
        else:
            logger.error("ChIP-seq only supported for bowtie2 and bwa.")
            sys.exit(-1)
        unique_bam = filterer(dd.get_work_bam(data), data)
        data["work_bam"] = unique_bam
    else:
        logger.info("Warning: When BAM file is given as input, bcbio skips multimappers removal."
                    "If BAM is not cleaned for peak calling, can result in downstream errors.")
    return [[data]]
