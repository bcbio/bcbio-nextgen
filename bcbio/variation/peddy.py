"""
correspondence checking with peddy (https://github.com/brentp/peddy)
"""
import collections
import os
import shutil
import sys

import toolz as tz

from collections import defaultdict
from bcbio.distributed.transaction import tx_tmpdir
from bcbio import utils
from bcbio.cwl import cwlutils
from bcbio.utils import safe_makedir, file_exists
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.log import logger
from bcbio.variation import vcfutils
from bcbio.variation.population import create_ped_file
from bcbio.qc import variant
from bcbio.provenance import do


PEDDY_OUT_EXTENSIONS = [".background_pca.json", ".het_check.csv", ".pca_check.png",
                        ".ped_check.png", ".ped_check.rel-difference.csv",
                        ".ped_check.csv", ".peddy.ped", ".sex_check.csv", ".ped_check.png",
                        ".html"]

def run_peddy_parallel(samples, parallel_fn):
    batch_samples = get_samples_by_batch(samples)
    to_run = batch_samples.values()
    samples = parallel_fn("run_peddy", [[x] for x in to_run])
    return [[utils.to_single_data(x)] for x in samples]

def run_qc(_, data, out_dir):
    """Run quality control in QC environment on a single sample.

    Enables peddy integration with CWL runs.
    """
    if cwlutils.is_cwl_run(data):
        qc_data = run_peddy([data], out_dir)
        if tz.get_in(["summary", "qc", "peddy"], qc_data):
            return tz.get_in(["summary", "qc", "peddy"], qc_data)

def is_human(data):
    return (tz.get_in(["genome_resources", "aliases", "human"], data, False) or
            dd.get_genome_build(data) in ["hg19", "GRCh37", "hg38"])

def run_peddy(samples, out_dir=None):
    vcf_file = None
    for d in samples:
        vcinfo = variant.get_active_vcinfo(d, use_ensemble=False)
        if vcinfo and vcinfo.get("vrn_file") and utils.file_exists(vcinfo["vrn_file"]):
            if vcinfo["vrn_file"] and dd.get_sample_name(d) in vcfutils.get_samples(vcinfo["vrn_file"]):
                vcf_file = vcinfo["vrn_file"]
                break
    data = samples[0]
    peddy = config_utils.get_program("peddy", data) if config_utils.program_installed("peddy", data) else None
    if not peddy or not vcf_file or not is_human(data):
        logger.info("peddy is not installed, not human or sample VCFs don't match, skipping correspondence checking "
                    "for %s." % vcf_file)
        return samples
    batch = dd.get_batch(data) or dd.get_sample_name(data)
    if out_dir:
        peddy_dir = safe_makedir(out_dir)
    else:
        peddy_dir = safe_makedir(os.path.join(dd.get_work_dir(data), "qc", batch, "peddy"))
    ped_file = create_ped_file(samples, vcf_file, out_dir=out_dir)
    peddy_prefix = os.path.join(peddy_dir, batch)
    peddy_report = peddy_prefix + ".html"
    peddyfiles = expected_peddy_files(peddy_report, batch)
    if file_exists(peddy_report):
        return dd.set_in_samples(samples, dd.set_summary_qc, peddyfiles)
    if file_exists(peddy_prefix + "-failed.log"):
        return samples
    num_cores = dd.get_num_cores(data)

    with tx_tmpdir(data) as tx_dir:
        peddy_prefix_tx = os.path.join(tx_dir, os.path.basename(peddy_prefix))
        # Redirects stderr because incredibly noisy with no intervals found messages from cyvcf2
        stderr_log = os.path.join(tx_dir, "run-stderr.log")
        sites_str = "--sites hg38" if dd.get_genome_build(data) == "hg38" else ""
        cmd = ("{peddy} -p {num_cores} {sites_str} --plot --prefix {peddy_prefix_tx} "
               "{vcf_file} {ped_file} 2> {stderr_log}")
        message = "Running peddy on {vcf_file} against {ped_file}."
        try:
            do.run(cmd.format(**locals()), message.format(**locals()))
        except:
            to_show = collections.deque(maxlen=100)
            with open(stderr_log) as in_handle:
                for line in in_handle:
                    to_show.append(line)
            def allowed_errors(l):
                return ((l.find("IndexError") >= 0 and l.find("is out of bounds for axis") >= 0) or
                        (l.find("n_components=") >= 0 and l.find("must be between 1 and n_features=") >= 0))
            def all_line_errors(l):
                return (l.find("no intervals found for") >= 0)
            if any([allowed_errors(l) for l in to_show]) or all([all_line_errors(l) for l in to_show]):
                logger.info("Skipping peddy because no variants overlap with checks: %s" % batch)
                with open(peddy_prefix + "-failed.log", "w") as out_handle:
                    out_handle.write("peddy did not find overlaps with 1kg sites in VCF, skipping")
                return samples
            else:
                logger.warning("".join(to_show))
                raise
        for ext in PEDDY_OUT_EXTENSIONS:
            if os.path.exists(peddy_prefix_tx + ext):
                shutil.move(peddy_prefix_tx + ext, peddy_prefix + ext)
    return dd.set_in_samples(samples, dd.set_summary_qc, peddyfiles)

def get_samples_by_batch(samples):
    batch_samples = defaultdict(list)
    for data in dd.sample_data_iterator(samples):
        batch = dd.get_batch(data) or dd.get_sample_name(data)
        if isinstance(batch, list):
            batch = tuple(batch)
        batch_samples[batch].append(data)
    return batch_samples

def expected_peddy_files(peddy_report, batch):
    peddydir = os.path.dirname(peddy_report)
    secondary = [batch + x for x in PEDDY_OUT_EXTENSIONS]
    secondary = [os.path.join(peddydir, x) for x in secondary]
    secondary = [x for x in secondary if os.path.exists(x)]
    return {"peddy": {"base": peddy_report, "secondary": secondary}}
