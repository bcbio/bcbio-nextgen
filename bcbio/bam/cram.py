"""Handle conversions to/from CRAM reference based compression.

http://www.ebi.ac.uk/ena/about/cram_toolkit
"""
import os
import subprocess

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
    resources = config_utils.get_resources("cram", config)
    jvm_opts = " ".join(resources.get("jvm_opts", ["-Xmx750m", "-Xmx2g"]))
    cram_jar = config_utils.get_jar("cramtools",
                                    config_utils.get_program("cram", config, "dir"))
    samtools = config_utils.get_program("samtools", config)
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            orig_header = "%s-header.sam" % os.path.splitext(out_file)[0]
            header_cmd = "{samtools} view -H -o {orig_header} {in_file}"
            cmd = ("java {jvm_opts} -jar {cram_jar} cram --input-bam-file {in_file} "
                   " --reference-fasta-file {ref_file} --preserve-read-names "
                   " --capture-all-tags --lossy-quality-score-spec '*8' "
                   "| java {jvm_opts} -jar {cram_jar} bam --output-bam-format "
                   "  --reference-fasta-file {ref_file} "
                   "| {samtools} reheader {orig_header} - "
                   "> {tx_out_file}")
            logger.info("Quality binning with CRAM")
            subprocess.check_call(header_cmd.format(**locals()), shell=True)
            subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file
