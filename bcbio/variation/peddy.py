"""
correspondence checking with peddy (https://github.com/brentp/peddy)
"""
import collections
import os
import shutil

import toolz as tz

from collections import defaultdict
from bcbio.distributed.transaction import tx_tmpdir
from bcbio import utils
from bcbio.utils import safe_makedir, file_exists
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.log import logger
from bcbio.variation import vcfutils
from bcbio.variation.population import create_ped_file
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

def run_peddy(samples):
    vcf_file = None
    for d in samples:
        if dd.get_vrn_file(d) and dd.get_sample_name(d) in vcfutils.get_samples(dd.get_vrn_file(d)):
            vcf_file = dd.get_vrn_file(d)
            break
    data = samples[0]
    peddy = config_utils.get_program("peddy", data) if config_utils.program_installed("peddy", data) else None
    is_human = tz.get_in(["genome_resources", "aliases", "human"], data, False)
    if not peddy or not vcf_file or not is_human:
        logger.info("peddy is not installed, not human or sample VCFs don't match, skipping correspondence checking "
                    "for %s." % vcf_file)
        return samples
    ped_file = create_ped_file(samples, vcf_file)
    batch = dd.get_batch(data) or dd.get_sample_name(data)
    peddy_dir = safe_makedir(os.path.join(dd.get_work_dir(data), "qc", batch, "peddy"))
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
        cmd = "{peddy} -p {num_cores} --plot --prefix {peddy_prefix_tx} {vcf_file} {ped_file} 2> {stderr_log}"
        message = "Running peddy on {vcf_file} against {ped_file}."
        try:
            do.run(cmd.format(**locals()), message.format(**locals()))
        except:
            to_show = collections.deque(maxlen=100)
            with open(stderr_log) as in_handle:
                for line in in_handle:
                    to_show.append(line)
            if to_show[-1].find("IndexError: index 0 is out of bounds for axis 0 with size 0") >= 0:
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
    return {"peddy": {"base": peddy_report, "secondary": secondary}}
