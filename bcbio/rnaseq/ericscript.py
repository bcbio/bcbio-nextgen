import os
import functools
from collections import namedtuple
from bcbio.distributed.transaction import file_transaction
import bcbio.pipeline.datadict as dd
from bcbio.provenance import do
from bcbio import utils


def run(config):
    ericscript = EricScriptConfig(config)

    with file_transaction(config, ericscript.output_dir) as tx_output_dir:
        input_data = get_input_data(tx_output_dir, config)
        if input_data.convert_cmd:
            msg = 'Convert disambiguated bam reads to fastq'
            do.run(input_data.convert_cmd, msg, env=ericscript.env)
        cmd = ericscript.get_run_command(tx_output_dir, input_data.files)
        do.run(cmd, ericscript.info_message, env=ericscript.env)


InputData = namedtuple('InputData', ['files', 'convert_cmd'])


def get_input_data(tx_output_dir, config):
    if dd.get_disambiguate(config):
        work_bam = dd.get_work_bam(config)
        fq_files = _get_fq_fnames(tx_output_dir)
        convert_cmd = [
            'bamToFastq',
            '-i', work_bam,
            '-fq', fq_files[0],
            '-fq2', fq_files[1]
        ]
    else:
        fq_files = dd.get_input_sequence_files(config)
        convert_cmd = None
    return InputData(fq_files, convert_cmd)


def _get_fq_fnames(location):
    join = functools.partial(os.path.join, location)
    return map(join, ('input1.fq', 'input2.fq'))


class EricScriptConfig(object):
    _OUTPUT_DIR_NAME = "ericscript"
    _MESSAGE = "Detect gene fusions with EricScript"

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
            '--remove',
        ] + list(input_files)

    @property
    def info_message(self):
        return self._MESSAGE

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
        return dd.get_ericscript_outdir(config)

    def _get_ericscript_db(self, config):
        return '/data/ericscript/ericscript_db'
