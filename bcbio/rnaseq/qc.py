"""Run Broad's RNA-SeqQC tool and handle reporting of useful summary metrics.
"""

import csv
import os

from bcbio import bam
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.utils import safe_makedir, file_exists
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


def sample_summary(bam_file, data, out_dir):
    """Run RNA-SeQC on a single RNAseq sample, writing to specified output directory.
    """
    metrics_file = os.path.join(out_dir, "metrics.tsv")
    if not file_exists(metrics_file):
        config = data["config"]
        ref_file = data["sam_ref"]
        genome_dir = os.path.dirname(os.path.dirname(ref_file))
        gtf_file = config_utils.get_transcript_gtf(genome_dir)
        rna_file = config_utils.get_rRNA_sequence(genome_dir)
        sample_file = os.path.join(safe_makedir(out_dir), "sample_file.txt")
        _write_sample_id_file(data, bam_file, sample_file)
        runner = rnaseqc_runner_from_config(config)
        broad_runner = runner_from_config(config)
        picard_index(broad_runner, bam_file)
        single_end = bam.is_paired(bam_file)
        runner.run(sample_file, ref_file, rna_file, gtf_file, out_dir, single_end)
    return _parse_rnaseqc_metrics(metrics_file, data["name"][-1])


def _write_sample_id_file(data, bam_file, out_file):
    HEADER = "\t".join(["Sample ID", "Bam File", "Notes"]) + "\n"
    sample_ids = ["\t".join([data["rgnames"]["pu"], bam_file, data["description"]])]
    with open(out_file, "w") as out_handle:
        out_handle.write(HEADER)
        for sample_id in sample_ids:
            out_handle.write(sample_id + "\n")
    return out_file

# ## Parsing

def _parse_rnaseqc_metrics(metrics_file, sample_name):
    """Parse RNA-SeQC tab delimited metrics file.
    """
    out = {}
    want = set(["Genes Detected", "Transcripts Detected",
                "Mean Per Base Cov.", "Estimated Library Size", "Fragment Length Mean",
                "Exonic Rate", "Intergenic Rate", "Intronic Rate",
                "Mapped", "Mapping Rate", "Duplication Rate of Mapped",
                "rRNA", "rRNA rate"])
    with open(metrics_file) as in_handle:
        reader = csv.reader(in_handle, dialect="excel-tab")
        header = reader.next()
        for metrics in reader:
            if metrics[1] == sample_name:
                for name, val in zip(header, metrics):
                    if name in want:
                        out[name] = val
    return out
