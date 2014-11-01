"""Identify structural variants using association testing with WHAM.

https://github.com/jewmanchue/wham
"""
import os
import subprocess
import sys

import toolz as tz
try:
    import vcf
except ImportError:
    vcf = None

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.variation import vcfutils
from bcbio.provenance import do

def run(items, background=None):
    """Detect copy number variations from batched set of samples using CNVkit.
    """
    if not background: background = []
    paired = vcfutils.get_paired_bams([x["align_bam"] for x in items], items)
    if paired:
        inputs = [paired.tumor_data]
        background_bams = [paired.normal_bam]
        background_names = [paired.normal_name]
    else:
        inputs = [items]
        background_bams = [x["align_bam"] for x in background]
        background_names = [dd.get_sample_name(x) for x in background]
    orig_vcf_file = _run_wham(inputs, background_bams)
    wclass_vcf_file = _add_wham_classification(orig_vcf_file, inputs)
    vcf_file = _fix_vcf(wclass_vcf_file, inputs, background_names)
    out = []
    for data in items:
        if "sv" not in data:
            data["sv"] = []
        data["sv"].append({"variantcaller": "wham", "vrn_file": vcf_file})
        out.append(data)
    return out

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "wham"))

def _run_wham(inputs, background_bams):
    """Run WHAM on a defined set of inputs and targets.
    """
    out_file = os.path.join(_sv_workdir(inputs[0]), "%s-wham.vcf" % dd.get_sample_name(inputs[0]))
    input_bams = [x["align_bam"] for x in inputs]
    if not utils.file_exists(out_file):
        with file_transaction(inputs[0], out_file) as tx_out_file:
            cores = dd.get_cores(inputs[0])
            background = "-b %s" % ",".join(background_bams) if background_bams else ""
            target_bams = ",".join(x["align_bam"] for x in inputs)
            target_bed = tz.get_in(["config", "algorithm", "variant_regions"], inputs[0])
            target_str = "-e %s" % target_bed if target_bed else ""
            cmd = ("WHAM-BAM -x {cores} -t {target_bams} {background} {target_str} "
                   "| sed 's/Numper/Number/' > {tx_out_file}")
            do.run(cmd.format(**locals()), "Run WHAM")
    return out_file

def _add_wham_classification(in_file, items):
    """Run WHAM classifier to assign a structural variant type to each call.
    """
    wham_sharedir = os.path.normpath(os.path.join(os.path.dirname(subprocess.check_output(["which", "WHAM-BAM"])),
                                                  os.pardir, "share", "wham"))
    out_file = "%s-class%s" % utils.splitext_plus(in_file)
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            this_python = sys.executable
            cmd = ("{this_python} {wham_sharedir}/classify_WHAM_vcf.py {in_file} "
                   "{wham_sharedir}/WHAM_training_data.txt > {tx_out_file}")
            do.run(cmd.format(**locals()), "Classify WHAM calls")
    return out_file

def _fix_vcf(orig_file, items, background_names):
    """Convert sample names and limit to larger structural calls.
    """
    out_file = "%s-clean%s" % utils.splitext_plus(orig_file)
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            reader = vcf.Reader(filename=orig_file)
            samples = [dd.get_sample_name(x) for x in items] + background_names
            assert len(samples) == len(reader.samples)
            reader.samples = samples
            writer = vcf.Writer(open(tx_out_file, "w"), template=reader)
            for rec in reader:
                if _is_sv(rec):
                    writer.write_record(rec)
    return out_file

def _is_sv(rec):
    """Pass along longer indels and structural variants.
    """
    size = max([len(x) for x in [rec.REF] + rec.ALT])
    is_sv = rec.is_sv or any(x for x in rec.ALT if str(x).startswith("<"))
    return size > 20 or is_sv
