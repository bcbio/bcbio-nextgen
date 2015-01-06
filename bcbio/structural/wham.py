"""Identify structural variants using association testing with WHAM.

https://github.com/jewmanchue/wham
"""
import csv
import os
import subprocess
import sys

import numpy as np
import toolz as tz
try:
    import vcf
except ImportError:
    vcf = None

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.structural import ensemble, shared
from bcbio.variation import bedutils, vcfutils
from bcbio.provenance import do

def run(items, background=None):
    """Detect copy number variations from batched set of samples using WHAM.
    """
    if not background: background = []
    paired = vcfutils.get_paired_bams([x["align_bam"] for x in items], items)
    if paired:
        inputs = [paired.tumor_data]
        background_bams = [paired.normal_bam]
        background_names = [paired.normal_name]
    else:
        assert not background
        inputs, background = shared.find_case_control(items)
        background_bams = [x["align_bam"] for x in background]
        background_names = [dd.get_sample_name(x) for x in background]
    orig_vcf_file = _run_wham(inputs, background_bams)
    wclass_vcf_file = _add_wham_classification(orig_vcf_file, inputs)
    vcf_file = _fix_vcf(wclass_vcf_file, inputs, background_names)
    bed_file = _convert_to_bed(vcf_file, inputs, use_lrt=len(background_bams) > 0)
    out = []
    for data in items:
        if "sv" not in data:
            data["sv"] = []
        data["sv"].append({"variantcaller": "wham",
                           "vrn_file": _subset_to_sample(bed_file, data),
                           "vcf_file": vcf_file})
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

def _subset_to_sample(bed_file, data):
    """Convert the global BED file into sample specific calls.
    """
    name = dd.get_sample_name(data)
    base, ext = os.path.splitext(bed_file)
    out_file = "%s-%s%s" % (base, name, ext)
    if not utils.file_uptodate(out_file, bed_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(bed_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    for line in in_handle:
                        sample_line = _check_bed_call(line, name)
                        if sample_line:
                            out_handle.write(sample_line)
    if utils.file_exists(out_file):
        bedprep_dir = utils.safe_makedir(os.path.join(os.path.dirname(out_file), "bedprep"))
        return bedutils.clean_file(out_file, data, bedprep_dir=bedprep_dir)
    else:
        return out_file

def _check_bed_call(line, name):
    """Check if a BED file line is called
    """
    chrom, start, end, stype_str, samples = line.split("\t")[:5]
    samples = set(samples.strip().split(";"))
    start, end = int(start), int(end)
    if name in samples:
        if end < start:
            fstart, fend = end, start
        else:
            fstart, fend = start, end
        if fend - fstart < ensemble.MAX_SVSIZE:
            return "%s\t%s\t%s\t%s\n" % (chrom, fstart, fend, stype_str)

def _convert_to_bed(vcf_file, inputs, use_lrt=False):
    """Convert WHAM output file into BED format for ensemble calling.

    Only outputs passing variants that have break end support and
    are above the mean of the likelihood ratio test
    score, if use_lrt is True.
    """
    buffer_size = 25  # bp around break ends
    out_file = "%s.bed" % utils.splitext_plus(vcf_file)[0]
    if not utils.file_uptodate(out_file, vcf_file):
        lrt_thresh = _calc_lrt_thresh(vcf_file) if use_lrt else 0.0
        with file_transaction(inputs[0], out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                writer = csv.writer(out_handle, dialect="excel-tab")
                for rec in _wham_pass_iter(vcf_file):
                    start = max(rec.start - buffer_size, 0)
                    _, end, _ = rec.INFO["BE"]
                    samples = [g.sample for g in rec.samples if g.gt_type > 0]
                    end = int(end) + buffer_size
                    if not use_lrt or rec.INFO.get("LRT", 0.0) > lrt_thresh:
                        writer.writerow([rec.CHROM, start, end, "%s_wham" % rec.INFO["WC"],
                                         ";".join(samples)])
    return out_file

def _calc_lrt_thresh(vcf_file):
    """Calculate threshold to use for including samples based on likelihood ratio test.

    For tumor/normal or case/control samples, this provides a threshold to include
    samples with higher likelihood in the cases. Calculates a simple mean
    threshold to use for inclusion.
    """
    lrts = [rec.INFO.get("LRT") for rec in _wham_pass_iter(vcf_file)]
    lrts = [x for x in lrts if x is not None]
    return np.mean(lrts)

def _wham_pass_iter(vcf_file):
    """Iterator over WHAM records that have breakend support and a called genotype.
    """
    reader = vcf.Reader(filename=vcf_file)
    for rec in reader:
        if rec.INFO["BE"][0] not in [".", None]:
            other_chrom, end, count = rec.INFO["BE"]
            if int(end) > int(rec.start) and other_chrom == rec.CHROM:
                samples = [g.sample for g in rec.samples if g.gt_type > 0]
                if len(samples) > 0:
                    yield rec

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
    return size > 10 or is_sv
