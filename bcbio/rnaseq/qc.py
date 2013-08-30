import os

from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.utils import safe_makedir
from bcbio.pipeline.qcsummary import is_paired
from bcbio.broad import runner_from_config
from bcbio.broad.picardrun import picard_index


class RNASeQCRunner(object):
    """
    Runs the Broad's RNA-SeQC tool:
    https://confluence.broadinstitute.org/display/CGATools/RNA-SeQC

    """
    def __init__(self, rnaseqc_path, bwa_path=None, jvm_opts=None):
        self._jvm_opts = " ".join(jvm_opts) if jvm_opts else "-Xms2g -Xmx4g"
        self._bwa_path = bwa_path if bwa_path else "bwa"
        self._rnaseqc_path = rnaseqc_path
        self._base_cmd = ("java -jar {jvm_opts} {rnaseqc_path} -n 1000 -s "
                          "{sample_file} -t {gtf_file} "
                          "-r {ref_file} -o {out_dir} -BWArRNA {rna_file} "
                          "-bwa {bwa_path} -ttype 2")

    def run(self, sample_file, ref_file, rna_file, gtf_file, out_dir,
            single_end=False):
        if single_end:
            self._base_cmd += " -singleEnd"
        cmd = self._base_cmd.format(rnaseqc_path=self._rnaseqc_path,
                                    bwa_path=self._bwa_path,
                                    jvm_opts=self._jvm_opts, **locals())
        do.run(cmd, "RNASeqQC on %s." % sample_file, None)


def rnaseqc_runner_from_config(config):
    """
    get a runner for Broad's RNA-SeQC tool using a bcbio-nextgen config dict to
    configure it
    """
    resources = config_utils.get_resources("rnaseqc", config)
    jvm_opts = resources.get("jvm_opts", ["-Xms750m", "-Xmx2g"])
    bwa_path = config_utils.get_program("bwa", config)
    rnaseqc_dir = config_utils.get_program("rnaseqc", config, "dir")
    rnaseqc_path = config_utils.get_jar("RNA-SeQC", rnaseqc_dir)
    return RNASeQCRunner(rnaseqc_path, bwa_path, jvm_opts)


def sample_summary(samples):
    sample_config = samples[0]
    config = sample_config[0]["config"]
    work_dir = sample_config[0]["dirs"]["work"]
    ref_file = sample_config[0]["sam_ref"]
    genome_dir = os.path.dirname(os.path.dirname(ref_file))
    gtf_file = config_utils.get_transcript_gtf(genome_dir)
    rna_file = config_utils.get_rRNA_sequence(genome_dir)

    out_dir = safe_makedir(os.path.join(work_dir, "qc", "rnaseqc"))
    sample_file = os.path.join(out_dir, "sample_file.txt")
    _write_sample_id_file(samples, sample_file)
    _index_samples(samples)
    runner = rnaseqc_runner_from_config(config)
    single_end = is_paired(sample_config[0]["work_bam"])
    runner.run(sample_file, ref_file, rna_file, gtf_file, out_dir, single_end)

    return samples


def _write_sample_id_file(samples, out_file):
    HEADER = "\t".join(["Sample ID", "Bam File", "Notes"]) + "\n"
    sample_ids = _extract_sample_ids(samples)
    with open(out_file, "w") as out_handle:
        out_handle.write(HEADER)
        for sample_id in sample_ids:
            out_handle.write(sample_id)
    return out_file


def _index_samples(samples):
    for data in samples:
        runner = runner_from_config(data[0]["config"])
        picard_index(runner, data[0]["work_bam"])


def _extract_sample_ids(samples):
    sample_ids = []
    for data in samples:
        names = data[0]["rgnames"]
        description = data[0]["description"]
        sample_ids.append("\t".join([names["pu"],
                                     data[0]["work_bam"], description]) + "\n")
    return sample_ids
