"""Runs gene fusion caller with EricScript.
Install EricScript via `bcbio upgrade --toolplus ericscript`,
or manually add the path to conda environment where it can be found
to the system config.
Reference data can be installed via `bcbio upgrade --datatarget ericscript`.
Alternatively, you can add path to the database into the system config.

EricScript requires bwa index to be built for its reference data.

To run gene fusion detection on disambiguated reads, we convert the .bam file
which was output by disambiguate to fastq files.

"""
import glob
import os

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.pipeline.fastq import convert_bam_to_fastq
from bcbio.provenance import do
from bcbio.pipeline import config_utils


def run(config):
    input_files = prepare_input_data(config)
    run_ericscript(config, input_files)
    return config


def prepare_input_data(config):
    """ In case of disambiguation, we want to run fusion calling on
    the disambiguated reads, which are in the work_bam file.
    As EricScript accepts 2 fastq files as input, we need to convert
    the .bam to 2 .fq files.
    """

    if not dd.get_disambiguate(config):
        return dd.get_input_sequence_files(config)

    work_bam = dd.get_work_bam(config)
    logger.info("Converting disambiguated reads to fastq...")
    fq_files = convert_bam_to_fastq(
        work_bam, dd.get_work_dir(config), None, None, config
    )
    return fq_files

def run_ericscript(data, input_files):
    db_location = dd.get_ericscript_db(data, None)
    if not db_location:
        logger.info("Skipping ericscript because ericscript database not found.")
    work_dir = dd.get_work_dir(data)
    sample_name = dd.get_sample_name(data)
    out_dir = os.path.join(work_dir, "ericscript", sample_name)
    ericscript = config_utils.get_program("ericscript.pl", data)
    ericscript_path = os.path.dirname(os.path.realpath(ericscript))
    samtools_path = os.path.join(ericscript_path, "..", "..", "bin")
    pathprepend = "export PATH=%s:$PATH; " % samtools_path
    files = " ".join(input_files)
    num_cores = dd.get_num_cores(data)
    cmd = ("{pathprepend} {ericscript} -db {db_location} -name {sample_name} -o {tx_out_dir} "
           "--nthreads {num_cores} {files}")
    message = "Running ericscript on %s using %s." %(files, db_location)
    with file_transaction(out_dir) as tx_out_dir:
        do.run(cmd.format(**locals()), message)
