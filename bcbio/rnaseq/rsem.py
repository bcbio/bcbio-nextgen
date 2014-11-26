import subprocess
from bcbio.distributed.transaction import tx_tmpdir, file_transaction
from bcbio.utils import chdir, which, safe_makedir
from bcbio import bam
from bcbio.log import logger
from bcbio.provenance import do

def prepare_rsem_reference(gtf, multifasta, build):
    """
    gtf: path to GTF file (must have gene_id and transcript_id)
    multifasta: path to multifasta file
    build: name of organism build (e.g. hg19)
    """
    if not which("rsem-prepare-reference"):
        logger.info("Skipping prepping RSEM reference because rsem-prepare-reference could "
                    "not be found.")
        return None

    cmd = "rsem-prepare-reference --gtf {gtf} {multifasta} {build}"
    with tx_tmpdir(remove=False) as rsem_genome_dir:
        with chdir(rsem_genome_dir):
            message = "Preparing rsem reference from %s" % gtf
            do.run(cmd.format(**locals()), message)
    return rsem_genome_dir

def rsem_calculate_expression(bam_file, rsem_genome_dir, samplename, build,
                              out_dir, cores=1):
    """
    works only in unstranded mode for now (--forward-prob 0.5)
    """
    if not which("rsem-calculate-expression"):
        logger.info("Skipping RSEM because rsem-calculate-expression could "
                    "not be found.")
        return None

    sentinel_file = os.path.join(out_dir, samplename + "Test.genes.results")
    if file_exists(sentinel_file):
        return out_dir

    paired_flag = "--paired" if bam.is_paired(bam_file) else ""
    core_flag = "-p {cores}".format(cores=cores)
    cmd = ("rsem-calculate-expression --bam {core_flag} {paired_flag} --no-bam-output "
           "--forward-prob 0.5 --estimate-rspd {bam_file} {rsem_genome_dir}/{build} "
           "{samplename}")
    message = "Calculating transcript expression of {bam_file} using RSEM."
    with file_transaction(out_dir) as tx_out_dir:
        safe_makedir(tx_out_dir)
        with chdir(tx_out_dir):
            do.run(cmd.format(**locals()), message.format(**locals()))
    return out_dir
