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
try:
    import pybedtools
except ImportError:
    pybedtools = None
import numpy
import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
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
    access_file = _create_access_file(dd.get_ref_file(items[0]), work_dir, items[0])
    if len(items + background) == 1:
        ckout = _run_cnvkit_single(items[0], access_file, work_dir)
    elif vcfutils.get_paired_phenotype(items[0]):
        ckout = _run_cnvkit_cancer(items, background, access_file, work_dir)
    else:
        ckout = _run_cnvkit_population(items, background, access_file, work_dir)
    ckout = _add_seg_to_output(ckout, items)
    ckout["variantcaller"] = "cnvkit"
    out = []
    for data in items:
        if "sv" not in data:
            data["sv"] = []
        data["sv"].append(ckout)
        out.append(data)
    return out

def _run_cnvkit_single(data, access_file, work_dir):
    """Process a single input file with a uniform background.
    """
    test_bams = [data["align_bam"]]
    background_bams = []
    return _run_cnvkit_shared(data, test_bams, background_bams, access_file, work_dir)

def _run_cnvkit_cancer(items, background, access_file, work_dir):
    """Run CNVkit on a tumor/normal pair.
    """
    paired = vcfutils.get_paired_bams([x["align_bam"] for x in items], items)
    return _run_cnvkit_shared(items[0], [paired.tumor_bam], [paired.normal_bam],
                              access_file, work_dir, background_name=paired.normal_name)

def _run_cnvkit_population(items, background, access_file, work_dir):
    """Run CNVkit on a population of samples.
    """
    raise NotImplementedError

def _run_cnvkit_shared(data, test_bams, background_bams, access_file, work_dir,
                       background_name=None):
    """Shared functionality to run CNVkit.
    """
    ref_file = dd.get_ref_file(data)
    raw_work_dir = os.path.join(work_dir, "raw")
    out_base = os.path.splitext(os.path.basename(test_bams[0]))[0]
    background_cnn = "%s_background.cnn" % (background_name if background_name else "flat")
    if not utils.file_exists(os.path.join(raw_work_dir, "%s.cnr" % out_base)):
        with tx_tmpdir(data, work_dir) as tx_work_dir:
            target_bed = tz.get_in(["config", "algorithm", "variant_regions"], data)
            cmd = ["batch"] + test_bams + ["-n"] + background_bams + ["-f", ref_file] + \
                  ["--targets", target_bed, "--access", access_file,
                   "-d", raw_work_dir, "--split",
                   "-p", str(tz.get_in(["config", "algorithm", "num_cores"], data, 1)),
                   "--output-reference", os.path.join(raw_work_dir, background_cnn)]
            at_avg, at_min = _get_antitarget_size(access_file, target_bed)
            if at_avg:
                cmd += ["--antitarget-avg-size", str(at_avg), "--antitarget-min-size", str(at_min)]
            args = cnvlib_cmd.parse_args(cmd)
            args.func(args)
            shutil.move(tx_work_dir, raw_work_dir)
    return {"cnr": os.path.join(raw_work_dir, "%s.cnr" % out_base),
            "cns": os.path.join(raw_work_dir, "%s.cns" % out_base),
            "back_cnn": os.path.join(raw_work_dir, background_cnn)}

def _add_seg_to_output(out, items):
    """Export outputs to 'seg' format compatible with IGV and GenePattern.
    """
    out_file = "%s.seg" % os.path.splitext(out["cns"])[0]
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            cmd = ["export", "seg", "-o", tx_out_file, out["cns"]]
            args = cnvlib_cmd.parse_args(cmd)
            args.func(args)
    out["seg"] = out_file
    return out

def _get_antitarget_size(access_file, target_bed):
    """Retrieve anti-target size based on distance between target regions.

    Handles smaller anti-target regions like found in subset genomes and tests.
    https://groups.google.com/d/msg/biovalidation/0OdeMfQM1CA/S_mobiz3eJUJ
    """
    prev = (None, 0)
    sizes = []
    for region in pybedtools.BedTool(access_file).subtract(target_bed):
        prev_chrom, prev_end = prev
        if region.chrom == prev_chrom:
            sizes.append(region.start - prev_end)
        prev = (region.chrom, region.end)
    avg_size = numpy.median(sizes) if len(sizes) > 0 else 0
    if avg_size < 1:
        return 2000, 200
    elif avg_size < 10000.0:  # Default antitarget-min-size
        return int(avg_size / 2.0), int(avg_size / 10.0)
    else:
        return None, None

def _create_access_file(ref_file, out_dir, data):
    """Create genome access file for CNVlib to define available genomic regions.

    XXX Can move to installation/upgrade process if too slow here.
    """
    out_file = os.path.join(out_dir, "%s-access.bed" % os.path.splitext(os.path.basename(ref_file))[0])
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [os.path.join(os.path.dirname(sys.executable), "genome2access.py"),
                   ref_file, "-s", "10000", "-o", tx_out_file]
            do.run(cmd, "Create CNVkit access file")
    return out_file
