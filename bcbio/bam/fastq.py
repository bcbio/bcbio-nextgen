from bcbio.log import logger
from bcbio import utils
from Bio import SeqIO


@utils.memoize_outfile(stem="groom")
def groom(in_file, in_qual="fastq-sanger", out_dir=None, out_file=None):
    """
    Grooms a FASTQ file into sanger format, if it is not already in that
    format. Use fastq-illumina for Illumina 1.3-1.7 qualities and
    fastq-solexa for the original solexa qualities. When in doubt, your
    sequences are probably fastq-sanger.

    """
    if in_qual == "fastq-sanger":
        logger.info("%s is already in Sanger format." % (in_file))
        return out_file

    count = SeqIO.convert(in_file, in_qual, out_file, "fastq-sanger")
    logger.info("Converted %d reads in %s to %s." % (count, in_file, out_file))
    return out_file
