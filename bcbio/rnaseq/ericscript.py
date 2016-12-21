import os
from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.pipeline.fastq import convert_bam_to_fastq
from bcbio.provenance import do


def run(config):
    # TODO build bwa index
    input_files = prepare_input_data(config)
    run_ericscript(config, input_files)
    return config


def prepare_input_data(config):
    """ Incase of disambiguation, we want to run fusion calling on
    the disambiguated reads, which are in the work_bam file.
    As EricScript accepts 2 fastq files as input, we need to convert
    the .bam to 2 .fq files.
    """

    if not dd.get_disambiguate(config):
        return dd.get_input_sequence_files(config)

    work_bam = dd.get_work_bam(config)
    fq_files = convert_bam_to_fastq(
        work_bam, dd.get_work_dir(config), None, None, config
    )
    return fq_files


def run_ericscript(sample_config, input_files):
    es_config = EricScriptConfig(sample_config)
    with file_transaction(sample_config, es_config.output_dir) as tx_out_dir:
        cmd = es_config.get_run_command(tx_out_dir, input_files)
        do.run(cmd, es_config.info_message, env=es_config.env)


class EricScriptConfig(object):
    info_message = "Detect gene fusions with EricScript"
    _OUTPUT_DIR_NAME = "ericscript"

    def __init__(self, config):
        self._env = self._get_env(config)
        self._sample_name = self._get_sample_name(config)
        self._db_location = self._get_ericscript_db(config)
        self._output_dir = self._get_output_dir(config)

    def get_run_command(self, tx_output_dir, input_files):
        return [
            'ericscript.pl',
            '-db', self._db_location,
            '-name', self._sample_name,
            '-o', tx_output_dir,
        ] + list(input_files)

    @property
    def env(self):
        return self._env

    @property
    def output_dir(self):
        return self._output_dir

    def _get_sample_name(self, config):
        return dd.get_lane(config)

    def _get_env(self, config):
        return utils.get_ericscript_env(config)

    def _get_output_dir(self, config):
        return os.path.join(dd.get_work_dir(config), self._OUTPUT_DIR_NAME)

    def _get_ericscript_db(self, config):
        # TODO get from config
        return '/data/ericscript/ericscript_db'
