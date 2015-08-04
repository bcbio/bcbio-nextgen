"""Run Broad's RNA-SeqQC tool and handle reporting of useful summary metrics.
"""

import csv
import os
from random import shuffle
from itertools import ifilter
import shutil
import uuid
import tempfile
import pandas as pd

# Provide transition period to install via upgrade with conda
try:
    import statsmodels.formula.api as sm
except ImportError:
    sm = None

from bcbio import bam
from bcbio import utils
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.utils import safe_makedir, file_exists
from bcbio.distributed.transaction import file_transaction
import bcbio.pipeline.datadict as dd

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
                          "-r {ref_file} -o {out_dir} -ttype 2 ")

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
        with file_transaction(data, out_dir) as tx_out_dir:
            config = data["config"]
            ref_file = data["sam_ref"]
            genome_dir = os.path.dirname(os.path.dirname(ref_file))
            gtf_file = dd.get_gtf_file(data)
            sample_file = os.path.join(safe_makedir(tx_out_dir), "sample_file.txt")
            _write_sample_id_file(data, bam_file, sample_file)
            runner = rnaseqc_runner_from_config(config)
            rna_file = config_utils.get_rRNA_sequence(genome_dir)
            bam.index(bam_file, config)
            single_end = not bam.is_paired(bam_file)
            runner.run(sample_file, ref_file, rna_file, gtf_file, tx_out_dir, single_end)
            # we don't need this large directory for just the report
            shutil.rmtree(os.path.join(tx_out_dir, data["description"]))
    return _parse_rnaseqc_metrics(metrics_file, data["name"][-1])

def _write_sample_id_file(data, bam_file, out_file):
    HEADER = "\t".join(["Sample ID", "Bam File", "Notes"]) + "\n"
    sample_ids = ["\t".join([data["description"], bam_file, data["description"]])]
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
                "Mean Per Base Cov.", "Fragment Length Mean",
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


def starts_by_depth(bam_file, data, sample_size=10000000):
    """
    Return a set of x, y points where x is the number of reads sequenced and
    y is the number of unique start sites identified
    If sample size < total reads in a file the file will be downsampled.
    """
    config = dd.get_config(data)
    binsize = (sample_size / 100) + 1
    seen_starts = set()
    counted = 0
    num_reads = []
    starts = []
    buffer = []
    downsampled = bam.downsample(bam_file, data, sample_size)
    with bam.open_samfile(downsampled) as samfile:
        for read in samfile:
            if read.is_unmapped:
                continue
            counted += 1
            buffer.append(str(read.tid) + ":" + str(read.pos))
            if counted % binsize == 0:
                seen_starts.update(buffer)
                buffer = []
                num_reads.append(counted)
                starts.append(len(seen_starts))
        seen_starts.update(buffer)
        num_reads.append(counted)
        starts.append(len(seen_starts))
    return pd.DataFrame({"reads": num_reads, "starts": starts})


def estimate_library_complexity(df, algorithm="RNA-seq"):
    """
    estimate library complexity from the number of reads vs.
    number of unique start sites. returns "NA" if there are
    not enough data points to fit the line
    """
    DEFAULT_CUTOFFS = {"RNA-seq": (0.25, 0.40)}
    cutoffs = DEFAULT_CUTOFFS[algorithm]
    if len(df) < 5:
        return {"unique_starts_per_read": 'nan',
                "complexity": "NA"}
    model = sm.ols(formula="starts ~ reads", data=df)
    fitted = model.fit()
    slope = fitted.params["reads"]
    if slope <= cutoffs[0]:
        complexity = "LOW"
    elif slope <= cutoffs[1]:
        complexity = "MEDIUM"
    else:
        complexity = "HIGH"

    # for now don't return the complexity flag
    return {"Unique Starts Per Read": float(slope)}
    # return {"unique_start_per_read": float(slope),
    #         "complexity": complexity}
