"""Identify structural variants using association testing with WHAM.

https://github.com/jewmanchue/wham
"""
import csv
import os
import shutil
import subprocess
import sys

import numpy as np
try:
    import vcf
except ImportError:
    vcf = None

from bcbio import utils
from bcbio.distributed.multi import run_multicore, zeromq_aware_logging
from bcbio.distributed.transaction import file_transaction
from bcbio.heterogeneity import chromhacks
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import region
from bcbio.structural import ensemble, shared
from bcbio.variation import bamprep, bedutils, vcfutils
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
    orig_bedpe = _run_wham(inputs, background_bams)
    #vcf_file = _fix_vcf(wclass_vcf_file, inputs, background_names)
    out = []
    for data in inputs:
        if "sv" not in data:
            data["sv"] = []
        data["sv"].append({"variantcaller": "wham",
                           "vrn_file": _get_sample_bed(orig_bedpe, data, background_names),
                           "vrn_bedpe": orig_bedpe})
        out.append(data)
    return out

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "wham"))

def _run_wham(inputs, background_bams):
    """Run WHAM on a defined set of inputs and targets.
    """
    out_file = os.path.join(_sv_workdir(inputs[0]), "%s-wham.bedpe" % dd.get_sample_name(inputs[0]))
    if not utils.file_exists(out_file):
        with file_transaction(inputs[0], out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                coords = chromhacks.autosomal_or_x_coords(dd.get_ref_file(inputs[0]))
                parallel = {"type": "local", "cores": dd.get_cores(inputs[0]), "progs": ["wham"]}
                rs = run_multicore(_run_wham_coords,
                                   [(inputs, background_bams, coord, out_file)
                                    for coord in coords],
                                   inputs[0]["config"], parallel)
                rs = {coord: fname for (coord, fname) in rs}
                for coord in coords:
                    with open(rs[coord]) as in_handle:
                        shutil.copyfileobj(in_handle, out_handle)
    return out_file

@utils.map_wrap
@zeromq_aware_logging
def _run_wham_coords(inputs, background_bams, coords, final_file):
    """Run WHAM on a specific set of chromosome, start, end coordinates.
    """
    base, ext = os.path.splitext(final_file)
    out_file = "%s-%s%s" % (base, region.to_safestr(coords), ext)
    if not utils.file_exists(out_file):
        with file_transaction(inputs[0], out_file) as tx_out_file:
            cores = dd.get_cores(inputs[0])
            ref_file = dd.get_ref_file(inputs[0])
            all_bams = ",".join([x["align_bam"] for x in inputs] + background_bams)
            coord_str = bamprep.region_to_gatk(coords)
            opts = "-k -m 30"
            cmd = ("WHAM-GRAPHENING {opts} -x {cores} -a {ref_file} -f {all_bams} -r {coord_str} "
                   "> {tx_out_file}")
            do.run(cmd.format(**locals()), "Run WHAM: %s" % region.to_safestr(coords))
    return [[coords, out_file]]

def _get_sample_bed(bedpe_file, data, background_samples):
    """Convert WHAM BEDPE output into BED format, selecting calls for the sample of interest.
    """
    background_samples = set(background_samples)
    sample = dd.get_sample_name(data)
    out_file = "%s-%s.bed" % (utils.splitext_plus(bedpe_file)[0], sample)
    if not utils.file_uptodate(out_file, bedpe_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(bedpe_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    reader = csv.reader(in_handle, dialect="excel-tab")
                    for parts in reader:
                        chrom1 = parts[0]
                        chrom2 = parts[3]
                        if chrom1 == chrom2:
                            info = {k: v for k, v in [x.split("=") for x in parts[-1].split(";")]}
                            cur_samples = set(info.get("LID", "").split(",") + info.get("RID", "").split(","))
                            if sample in cur_samples and len(cur_samples & background_samples) == 0:
                                start, end = info["POS"].split(",")
                                out_handle.write("%s\t%s\t%s\t%s\n" % (chrom1, start, end,
                                                                       "%s_wham" % info["SVTYPE"].strip()))
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
