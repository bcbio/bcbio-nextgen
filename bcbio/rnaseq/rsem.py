import os

from bcbio import bam
from bcbio import utils
from bcbio.log import logger
from bcbio.distributed import transaction
from bcbio.provenance import do

CALCULATE_EXP = (
    "rsem-calculate-expression --bam {core_flag} {paired_flag} "
    "--no-bam-output --forward-prob 0.5 "
    "--estimate-rspd {bam_file} {rsem_genome_dir}/{build} {samplename}")
PREPARE_REFERENCE = "rsem-prepare-reference --gtf {gtf} {multifasta} {build}"


def prepare_rsem_reference(gtf, multifasta, build):
    """
    gtf: path to GTF file (must have gene_id and transcript_id)
    multifasta: path to multifasta file
    build: name of organism build (e.g. hg19)
    """
    if not utils.which("rsem-prepare-reference"):
        logger.info("Skipping prepping RSEM reference because "
                    "rsem-prepare-reference could not be found.")
        return None

    command = PREPARE_REFERENCE.format(gtf=gtf, multifasta=multifasta,
                                       build=build)
    with transaction.tx_tmpdir(remove=False) as rsem_genome_dir:
        with utils.chdir(rsem_genome_dir):
            message = "Preparing rsem reference from %s" % gtf
            do.run(command, message)
    return rsem_genome_dir


def rsem_calculate_expression(bam_file, rsem_genome_dir, samplename,
                              build, out_dir, cores=1):
    """
    works only in unstranded mode for now (--forward-prob 0.5)
    """
    if not utils.which("rsem-calculate-expression"):
        logger.info("Skipping RSEM because rsem-calculate-expression could "
                    "not be found.")
        return None

    sentinel_file = os.path.join(out_dir, samplename + "Test.genes.results")
    if utils.file_exists(sentinel_file):
        return out_dir

    paired_flag = "--paired" if bam.is_paired(bam_file) else ""
    core_flag = "-p {cores}".format(cores=cores)
    command = CALCULATE_EXP.format(
        core_flag=core_flag, paired_flag=paired_flag, bam_file=bam_file,
        rsem_genome_dir=rsem_genome_dir, build=build, samplename=samplename)
    message = "Calculating transcript expression of {bam_file} using RSEM."

    with transaction.file_transaction(out_dir) as tx_out_dir:
        utils.safe_makedir(tx_out_dir)
        with utils.chdir(tx_out_dir):
            do.run(command, message.format(bam_file=bam_file))
    return out_dir
