import os

from bcbio.pipeline import config_utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.utils import (safe_makedir, file_exists)
from bcbio.provenance import do
from bcbio import utils
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.ngsalign import postalign
from bcbio import bam
from bcbio import broad


def align(fastq_file, pair_file, ref_file, names, align_dir, data):
    assert data["analysis"].lower().startswith("wgbs-seq"), "No comparible alignment"
    config = data["config"]
    sample = dd.get_sample_name(data)
    out_prefix = os.path.join(align_dir, dd.get_lane(data))
    ref_file = dd.get_sam_ref(data)

    final_out = os.path.join(align_dir, "{0}.bam".format(sample))
    if file_exists(final_out):
        data = dd.set_work_bam(data, final_out)
        return data

    bsmap = config_utils.get_program("bsmap", config)
    fastq_files = " -a %s" % fastq_file
    num_cores = dd.get_num_cores(data)
    num_cores = "-p %d" % num_cores
    safe_makedir(align_dir)
    cmd = "{bsmap} {num_cores} -w 100 -v 0.07 -m 10 -x 300 -o {tx_out_bam} -d {ref_file} {fastq_files}"
    if pair_file:
        fastq_files = "-a %s -b %s" % (fastq_file, pair_file)
    if not final_out:
        with file_transaction(final_out) as tx_out_bam:
            run_message = "Running BSMAP aligner on %s and %s" % (fastq_file, ref_file)
            do.run(cmd.format(**locals()), run_message, None)
    data = dd.set_work_bam(data, final_out)
    return data


def _process_bam(bam_file, in_fastq, sample, reference, config):
    broad_runner = broad.runner_from_config(config)
    names = {'rg': in_fastq, 'library': 'WGBS_LIB', 'pl': 'Illumina', 'pu': 'R1', 'sm': in_fastq, 'sample': sample}
    out_fix_bam = broad_runner.run_fn("picard_fix_rgs", bam_file, names)
    order_bam = utils.append_stem(out_fix_bam, "_order")
    broad_runner.run_fn("picard_reorder", out_fix_bam, reference, order_bam)
    bam.index(order_bam, config)
    # order_bam = _set_quality(order_bam)
    # bam.index(order_bam, config)
    return order_bam
