"""
correspondence checking with peddy (https://github.com/brentp/peddy)
"""

import os
from collections import defaultdict
from bcbio.utils import safe_makedir, file_exists
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.log import logger
from bcbio.variation.population import create_ped_file
from bcbio.provenance import do
from bcbio.distributed.transaction import tx_tmpdir

PEDDY_OUT_EXTENSIONS = [".background_pca.json", ".het_check.csv", ".pca_check.png",
                        ".ped_check.png", ".ped_check.rel-difference.csv",
                        ".peddy.ped", ".sex_check.csv", ".ped_check.png"]

def run_peddy_parallel(samples, parallel_fn):
    batch_samples = get_samples_by_batch(samples)
    to_run = batch_samples.values()
    samples = parallel_fn("run_peddy", [to_run])
    return samples

def run_peddy(samples):
    data = samples[0]
    vcf_file = dd.get_vrn_file(data)
    ped_file = create_ped_file(samples, vcf_file)
    peddy = config_utils.get_program("peddy", data)
    peddy_dir = os.path.join(dd.get_work_dir(data), "peddy", dd.get_batch(data))
    safe_makedir(peddy_dir)
    peddy_prefix = os.path.join(peddy_dir, dd.get_batch(data))
    peddy_report = peddy_prefix + ".html"
    batch = dd.get_batch(data)
    peddyfiles = expected_peddy_files(peddy_report, batch)
    if not peddy:
        logger.info("peddy is not installed, skipping correspondence checking "
                    "for %s." % vcf_file)
        return samples
    if file_exists(peddy_report):
        return dd.set_in_samples(samples, dd.set_summary_qc, peddyfiles)

    cmd = "{peddy} --plot --prefix {peddy_prefix} {vcf_file} {ped_file}"
    message = "Running peddy on {vcf_file} against {ped_file}."
    do.run(cmd.format(**locals()), message.format(**locals()))
    return dd.set_in_samples(samples, dd.set_summary_qc, peddyfiles)

def get_samples_by_batch(samples):
    batch_samples = defaultdict(list)
    for data in dd.sample_data_iterator(samples):
        batch = dd.get_batch(data)
        batch_samples[batch].append(data)
    return batch_samples

def expected_peddy_files(peddy_report, batch):
    peddydir = os.path.dirname(peddy_report)
    secondary = [batch + x for x in PEDDY_OUT_EXTENSIONS]
    secondary = [os.path.join(peddydir, x) for x in secondary]
    return {"peddy": {"base": peddy_report, "secondary": secondary}}
