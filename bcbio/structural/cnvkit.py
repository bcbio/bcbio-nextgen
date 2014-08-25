"""Copy number detection with CNVkit with specific support for targeted sequencing.

http://cnvkit.readthedocs.org
"""
import os
import shutil
import sys

try:
    from cnvlib import commands as cnvlib_cmd
except ImportError:
    cnvlib_cmd = None

import toolz as tz
from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.variation import vcfutils
from bcbio.provenance import do

def run(items, background=None):
    """Detect copy number variations from batched set of samples using CNVkit.
    """
    if not cnvlib_cmd:
        raise ImportError("cnvkit not installed")
    if not background: background = []
    work_dir = utils.safe_makedir(os.path.join(items[0]["dirs"]["work"], "structural",
                                               tz.get_in(["rgnames", "sample"], items[0]),
                                               "cnvkit"))
    return _cnvkit_by_type(items, background, work_dir)

def _cnvkit_by_type(items, background, work_dir):
    """Dispatch to specific CNVkit functionality based on input type.
    """
    if len(items + background) == 1:
        return _run_cnvkit_single(items[0], work_dir)
    elif vcfutils.get_paired_phenotype(items[0]):
        return _run_cnvkit_cancer(items, background, work_dir)
    else:
        return _run_cnvkit_population(items, background, work_dir)

def _run_cnvkit_single(data, work_dir):
    """Process a single input file with a uniform background.
    """
    ref_file = dd.get_ref_file(data)
    access_file = _create_access_file(ref_file, work_dir)
    raw_work_dir = os.path.join(work_dir, "raw")
    if not utils.file_exists(os.path.join(raw_work_dir, "%s.cnr" %
                                          os.path.splitext(os.path.basename(data["align_bam"]))[0])):
        with utils.curdir_tmpdir(data, work_dir) as tx_work_dir:
            cmd = ["batch", data["align_bam"], "-n", "-f", ref_file,
                   "--targets", tz.get_in(["config", "algorithm", "variant_regions"], data),
                   "--access", access_file,
                   "-d", tx_work_dir, "--split", "-p", str(tz.get_in(["config", "algorithm", "num_cores"], data, 1)),
                   # XXX Currently set for testing -- need to generalize
                   "--target-avg-size", "500", "--antitarget-avg-size", "500"]
            args = cnvlib_cmd.parse_args(cmd)
            args.func(args)
            shutil.move(tx_work_dir, raw_work_dir)
    raise NotImplementedError

def _create_access_file(ref_file, out_dir):
    """Create genome access file for CNVlib to define available genomic regions.

    XXX Need to move to installation/upgrade process.
    """
    out_file = os.path.join(out_dir, "%s-access.bed" % os.path.splitext(os.path.basename(ref_file))[0])
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            cmd = [os.path.join(os.path.dirname(sys.executable), "genome2access.py"),
                   ref_file, "-s", "10000", "-o", tx_out_file]
            do.run(cmd, "Create CNVkit access file")
    return out_file

def _run_cnvkit_cancer(items, background, work_dir):
    raise NotImplementedError

def _run_cnvkit_population(items, background, work_dir):
    raise NotImplementedError
