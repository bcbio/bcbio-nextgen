import os

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.ngsalign import bwa
from bcbio.pipeline import datadict as dd
from bcbio.pipeline.fastq import convert_bam_to_fastq
from bcbio.provenance import do


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


def run_ericscript(sample_config, input_files):
    es_config = EricScriptConfig(sample_config)
    utils.safe_makedir(es_config.output_dir)

    with file_transaction(sample_config, es_config.sample_out_dir) as tx_out:
        cmd = es_config.get_run_command(tx_out, input_files)
        logger.info("Running EricScript:\n%s" % ' '.join(cmd))
        do.run(cmd, es_config.info_message, env=es_config.env)


def build_bwa_index_if_absent(es_config, sample_config):
    if not os.path.exists(es_config.reference_index):
        bwa.build_bwa_index(es_config.reference_fasta, sample_config)


class EricScriptConfig(object):
    info_message = 'Detect gene fusions with EricScript'
    EXECUTABLE = 'ericscript.pl'
    _OUTPUT_DIR_NAME = 'ericscript'
    _REF_INDEX = 'data/homo_sapiens/allseq.fa.bwt'
    _REF_FASTA = 'data/homo_sapiens/allseq.fa'

    def __init__(self, config):
        self._db_location = dd.get_ericscript_db(config)
        self._env_prefix = dd.get_ericscript_env(config)
        self._sample_name = dd.get_lane(config)
        self._work_dir = dd.get_work_dir(config)
        self._env = None
        self._output_dir = None
        self._sample_out_dir = None

    def get_run_command(self, tx_output_dir, input_files):
        logger.debug("Input data: %s" % ', '.join(input_files))
        return [
            self.EXECUTABLE,
            '-db', self._db_location,
            '-name', self._sample_name,
            '-o', tx_output_dir,
        ] + list(input_files)

    @property
    def env(self):
        if self._env is None:
            self._env = self._get_env()
        return self._env

    @property
    def output_dir(self):
        if self._output_dir is None:
            self._output_dir = self._get_output_dir()
        return self._output_dir

    @property
    def sample_out_dir(self):
        if self._sample_out_dir is None:
            self._sample_out_dir = os.path.join(
                self.output_dir, self._sample_name
            )
        return self._sample_out_dir

    @property
    def reference_index(self):
        return os.path.join(self._db_location, self._REF_INDEX)

    @property
    def reference_fasta(self):
        return os.path.join(self._db_location, self._REF_FASTA)

    def _get_env(self):
        return utils.get_ericscript_env(self._env_prefix)

    def _get_output_dir(self):
        return os.path.join(self._work_dir, self._OUTPUT_DIR_NAME)
