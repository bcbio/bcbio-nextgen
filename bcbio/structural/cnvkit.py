"""Copy number detection with CNVkit with specific support for targeted sequencing.

http://cnvkit.readthedocs.org
"""
import os
import shutil
import sys

try:
    import pybedtools
except ImportError:
    pybedtools = None
import numpy
import toolz as tz

from bcbio import install, utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.pipeline import datadict as dd
from bcbio.variation import vcfutils
from bcbio.provenance import do
from bcbio.structural import shared, theta

def run(items, background=None):
    """Detect copy number variations from batched set of samples using CNVkit.
    """
    if not background: background = []
    return _cnvkit_by_type(items, background)

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "cnvkit"))

def _cnvkit_by_type(items, background):
    """Dispatch to specific CNVkit functionality based on input type.
    """
    if len(items + background) == 1:
        return _run_cnvkit_single(items[0])
    elif vcfutils.get_paired_phenotype(items[0]):
        return _run_cnvkit_cancer(items, background)
    else:
        return _run_cnvkit_population(items, background)

def _associate_cnvkit_out(ckout, items):
    """Associate cnvkit output with individual items.
    """
    ckout = _add_seg_to_output(ckout, items[0])
    ckout["variantcaller"] = "cnvkit"
    out = []
    for data in items:
        ckout = _add_bed_to_output(ckout, data)
        if "sv" not in data:
            data["sv"] = []
        data["sv"].append(ckout)
        out.append(data)
    return out

def _run_cnvkit_single(data, access_file=None, background=None):
    """Process a single input file with BAM or uniform background.
    """
    work_dir = _sv_workdir(data)
    if not access_file:
        access_file = _create_access_file(dd.get_ref_file(data), work_dir, data)
    test_bams = [data["align_bam"]]
    if background:
        background_bams = [x["align_bam"] for x in background]
        background_name = os.path.splitext(os.path.basename(background_bams[0]))[0]
    else:
        background_bams = []
        background_name = None
    ckout = _run_cnvkit_shared(data, test_bams, background_bams, access_file, work_dir,
                               background_name=background_name)
    return _associate_cnvkit_out(ckout, [data])

def _run_cnvkit_cancer(items, background):
    """Run CNVkit on a tumor/normal pair.
    """
    paired = vcfutils.get_paired_bams([x["align_bam"] for x in items], items)
    work_dir = _sv_workdir(paired.tumor_data)
    access_file = _create_access_file(dd.get_ref_file(paired.tumor_data), work_dir, paired.tumor_data)
    ckout = _run_cnvkit_shared(paired.tumor_data, [paired.tumor_bam], [paired.normal_bam],
                               access_file, work_dir, background_name=paired.normal_name)
    ckout = theta.run(ckout, paired)
    tumor_data = _associate_cnvkit_out(ckout, [paired.tumor_data])
    normal_data = [x for x in items if dd.get_sample_name(x) != paired.tumor_name]
    return tumor_data + normal_data

def _run_cnvkit_population(items, background):
    """Run CNVkit on a population of samples.

    Tries to calculate background based on case/controls, otherwise uses
    a flat background for each sample and calls independently.
    """
    assert not background
    inputs, background = shared.find_case_control(items)
    access_file = _create_access_file(dd.get_ref_file(inputs[0]), _sv_workdir(inputs[0]), inputs[0])
    return [_run_cnvkit_single(data, access_file, background)[0] for data in inputs] + \
           [_run_cnvkit_single(data, access_file, inputs)[0] for data in background]

def _run_cnvkit_shared(data, test_bams, background_bams, access_file, work_dir,
                       background_name=None):
    """Shared functionality to run CNVkit.
    """
    ref_file = dd.get_ref_file(data)
    raw_work_dir = os.path.join(work_dir, "raw")
    out_base = os.path.splitext(os.path.basename(test_bams[0]))[0]
    background_cnn = "%s_background.cnn" % (background_name if background_name else "flat")
    if not utils.file_exists(os.path.join(raw_work_dir, "%s.cnr" % out_base)):
        if os.path.exists(raw_work_dir):
            shutil.rmtree(raw_work_dir)
        with tx_tmpdir(data, work_dir) as tx_work_dir:
            target_bed = tz.get_in(["config", "algorithm", "variant_regions"], data)
            cores = min(tz.get_in(["config", "algorithm", "num_cores"], data, 1),
                        len(test_bams) + len(background_bams))
            cmd = [os.path.join(os.path.dirname(sys.executable), "cnvkit.py"), "batch"] + \
                  test_bams + ["-n"] + background_bams + ["-f", ref_file] + \
                  ["--targets", target_bed, "--access", access_file,
                   "-d", tx_work_dir, "--split", "-p", str(cores),
                   "--output-reference", os.path.join(tx_work_dir, background_cnn)]
            at_avg, at_min, t_avg = _get_antitarget_size(access_file, target_bed)
            if at_avg:
                cmd += ["--antitarget-avg-size", str(at_avg), "--antitarget-min-size", str(at_min),
                        "--target-avg-size", str(t_avg)]
            local_sitelib = os.path.join(install.get_defaults().get("tooldir", "/usr/local"),
                                         "lib", "R", "site-library")
            cmd += ["--rlibpath", local_sitelib]
            do.run(cmd, "CNVkit batch")
            shutil.move(tx_work_dir, raw_work_dir)
    return {"cnr": os.path.join(raw_work_dir, "%s.cnr" % out_base),
            "cns": os.path.join(raw_work_dir, "%s.cns" % out_base),
            "back_cnn": os.path.join(raw_work_dir, background_cnn)}

def _add_seg_to_output(out, data):
    """Export outputs to 'seg' format compatible with IGV and GenePattern.
    """
    out_file = "%s.seg" % os.path.splitext(out["cns"])[0]
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [os.path.join(os.path.dirname(sys.executable), "cnvkit.py"), "export",
                   "seg", "-o", tx_out_file, out["cns"]]
            do.run(cmd, "CNVkit export seg")
    out["seg"] = out_file
    return out

def _add_bed_to_output(out, data):
    """Add FreeBayes cnvmap BED-like representation to the output.
    """
    out_file = "%s.bed" % os.path.splitext(out["cns"])[0]
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [os.path.join(os.path.dirname(sys.executable), "cnvkit.py"), "export",
                   "freebayes", "--name", dd.get_sample_name(data),
                   "--ploidy", str(dd.get_ploidy(data)),
                   "-o", tx_out_file, out["cns"]]
            gender = dd.get_gender(data)
            if gender:
                cmd += ["--gender", gender]
            do.run(cmd, "CNVkit export FreeBayes BED cnvmap")
    out["vrn_file"] = out_file
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
    if avg_size < 10000.0:  # Default antitarget-min-size
        return 1000, 75, 1000
    else:
        return None, None, None

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
