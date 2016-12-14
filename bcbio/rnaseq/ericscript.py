import os
from bcbio.distributed.transaction import file_transaction
import bcbio.pipeline.datadict as dd
from bcbio.provenance import do
from bcbio import utils


def run(config):
    env = utils.get_ericscript_env(config)
    db_location = _get_ericscript_db(config)
    fastq_files = _get_fastq_files(config)
    samplename = dd.get_lane(config)
    output_location = get_output_dir(config)

    with file_transaction(config, output_location) as tx_output_location:
        cmd = [
            'ericscript.pl',
            '-db', db_location,
            '-name', samplename,
            '-o', tx_output_location,
        ]
        cmd += fastq_files
        do.run(cmd, "Detect gene fusions with EricScript", env=env)


def get_output_dir(config):
    return os.path.join(dd.get_work_dir(config), "ericscript")


def _get_ericscript_db(config):
    return '/data/ericscript/ericscript_db'


def _get_fastq_files(config):
    return dd.get_input_sequence_files(config)
