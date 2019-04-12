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
    es_config = EricScriptConfig(data)
    if not es_config.has_ericscript_db():
        logger.info("Skipping ericscript because ericscript database not found.")
        return None
    db_location = es_config._db_location
    work_dir = dd.get_work_dir(data)
    sample_name = dd.get_sample_name(data)
    out_dir = os.path.join(work_dir, "ericscript", sample_name)
    ericscript = config_utils.get_program("ericscript.pl", data)
    ericscript_path = os.path.dirname(os.path.realpath(ericscript))
    samtools_path = os.path.join(ericscript_path, "..", "..", "bin")
    pathprepend = "export PATH=%s:\"$PATH\"; " % samtools_path
    files = " ".join(input_files)
    num_cores = dd.get_num_cores(data)
    cmd = ("{pathprepend} {ericscript} -db {db_location} -name {sample_name} -o {tx_out_dir} "
           "--nthreads {num_cores} {files}")
    message = "Running ericscript on %s using %s." %(files, db_location)
    with file_transaction(out_dir) as tx_out_dir:
        do.run(cmd.format(**locals()), message)

class EricScriptConfig(object):
    """This class which encapsulates access to the dat
    related to EricScript in the sample config dictionary.

    Public constants:
        info_message: text message passed as an argument to do.run
        EXECUTABLE: name of the EricScipt command

    Private constants:
        _OUTPUT_DIR_NAME: name of the dir created in working directory for
    ericscript output
    """
    info_message = 'Detect gene fusions with EricScript'
    EXECUTABLE = 'ericscript.pl'
    _OUTPUT_DIR_NAME = 'ericscript'
    _REF_INDEX = 'allseq.fa.bwt'
    _REF_FASTA = 'allseq.fa'

    def __init__(self, data):
        self._db_location = self._get_ericscript_db(data)
        self._sample_name = dd.get_lane(data)
        self._work_dir = dd.get_work_dir(data)
        self._env = None
        self._output_dir = None
        self._sample_out_dir = None

    def _get_ericscript_db(self, data):
        transcript_file = dd.get_gtf_file(data)
        if transcript_file and os.path.exists(transcript_file):
            transcript_dir = os.path.dirname(transcript_file)
            ericscript_dirs = glob.glob(os.path.join(transcript_dir, "ericscript", "ericscript_db*"))
            if ericscript_dirs:
                return sorted(ericscript_dirs)[-1]

    def has_ericscript_db(self):
        return self._db_location is not None

    def get_run_command(self, tx_output_dir, input_files):
        """Constructs a command to run EricScript via do.run function.

        :param tx_output_dir: A location where all EricScript output will be
        written during execution.
        :param input_files: an iterable with paths to 2 fastq files
        with input data.
        :return: list
        """
        logger.debug("Input data: %s" % ', '.join(input_files))
        cmd = [
            self.EXECUTABLE,
            '-db', self._db_location,
            '-name', self._sample_name,
            '-o', tx_output_dir,
        ] + list(input_files)
        return "export PATH=%s:%s:\"$PATH\"; %s;" % (self._get_samtools0_path(), self._get_ericscript_path(), " ".join(cmd))

    def _get_ericscript_path(self):
        """Retrieve PATH to the isolated eriscript anaconda environment.
        """
        es = utils.which(os.path.join(utils.get_bcbio_bin(), self.EXECUTABLE))
        return os.path.dirname(os.path.realpath(es))
    def _get_samtools0_path(self):
        """Retrieve PATH to the samtools version specific for eriscript.
        """
        samtools_path = os.path.realpath(os.path.join(self._get_ericscript_path(),"..", "..", "bin"))
        return samtools_path

    @property
    def output_dir(self):
        """Absolute path to permanent location in working directory
        where EricScript output will be stored.
        """
        if self._output_dir is None:
            self._output_dir = self._get_output_dir()
        return self._output_dir

    @property
    def sample_out_dir(self):
        """Absolute path to permanent location in working directory
        where EricScript output for the current sample will be stored.
        (a subdirectory of `output_dir`)
        """
        if self._sample_out_dir is None:
            self._sample_out_dir = os.path.join(
                self.output_dir, self._sample_name
            )
        return self._sample_out_dir

    @property
    def reference_index(self):
        """Absolute path to the BWA index for EricScript reference data."""
        if self._db_location:
            ref_indices = glob.glob(os.path.join(self._db_location, "*", self._REF_INDEX))
            if ref_indices:
                return ref_indices[0]

    @property
    def reference_fasta(self):
        """Absolute path to the fasta file with EricScript reference data."""
        if self._db_location:
            ref_files = glob.glob(os.path.join(self._db_location, "*", self._REF_FASTA))
            if ref_files:
                return ref_files[0]

    def _get_output_dir(self):
        return os.path.join(self._work_dir, self._OUTPUT_DIR_NAME)
