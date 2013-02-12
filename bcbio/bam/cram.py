"""Handle conversions to/from CRAM reference based compression.

http://www.ebi.ac.uk/ena/about/cram_toolkit
"""
import os
import sh
from bcbio.log import logger
from bcbio.pipeline import config_utils
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction


def illumina_qual_bin(in_file, ref_file, out_dir, config):
    """Uses CRAM to perform Illumina 8-bin approaches to existing BAM files.

    Bins quality scores according to Illumina scheme:

    http://www.illumina.com/Documents/products/whitepapers/whitepaper_datacompression.pdf

    Also fixes output header to remove extra run groups added by CRAM during conversion.
    """
    index_file = ref_file + ".fai"
    assert os.path.exists(index_file), "Could not find FASTA reference index: %s" % index_file
    out_file = os.path.join(out_dir, "%s-qualbin%s" % os.path.splitext(os.path.basename(in_file)))
    cram_jar = config_utils.get_jar("cramtools",
                                    config_utils.get_program("cram", config, "dir"))
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            orig_header = "%s-header.sam" % os.path.splitext(out_file)[0]
            java_cmd = sh.Command("java").bake("-jar", cram_jar, _long_sep=" ")
            to_cram = java_cmd.bake("cram", input_bam_file=in_file,
                                    reference_fasta_file=ref_file,
                                    preserve_read_names=True,
                                    capture_all_tags=True,
                                    lossy_quality_score_spec="'*8'")
            to_bam = java_cmd.bake("bam", output_bam_format=True,
                                   reference_fasta_file=ref_file)
            make_header = sh.samtools.view.bake(H=True, o=orig_header)
            make_header(in_file)
            reheader = sh.samtools.reheader.bake(orig_header, "-")
            logger.info("Quality binning with CRAM")
            reheader(to_bam(to_cram()), _out=tx_out_file)
    return out_file
