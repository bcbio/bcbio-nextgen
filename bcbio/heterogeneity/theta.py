"""Provide estimates of sample purity and subclonal copy number using THetA.

Identifying cellularity and subclonal populations within somatic calling using
tumor normal pairs.

https://github.com/raphael-group/THetA
"""
import os
import sys
import subprocess

import pybedtools
import pysam
import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.variation import vcfutils
from bcbio.provenance import do
from bcbio.structural import cnvkit, convert

def run(vrn_info, cnvs_by_name, somatic_info):
    """Run THetA analysis given output from CNV caller on a tumor/normal pair.
    """
    cmd = _get_cmd("RunTHeTA.py")
    if not cmd:
        logger.info("THetA scripts not found in current PATH. Skipping.")
    else:
        work_dir = _sv_workdir(somatic_info.tumor_data)
        assert "cnvkit" in cnvs_by_name, "THetA requires CNVkit calls"
        cnv_info = cnvkit.export_theta(cnvs_by_name["cnvkit"], somatic_info.tumor_data)
        cnv_info["theta_input"] = subset_by_supported(cnv_info["theta_input"], _theta_to_coords,
                                                      cnvs_by_name, work_dir, somatic_info.tumor_data)
        return _run_theta(cnv_info, somatic_info.tumor_data, work_dir, run_n3=False)

def _theta_to_coords(line):
    _, chrom, start, end = line.split()[:4]
    return (chrom, start, end)

def subset_by_supported(input_file, get_coords, calls_by_name, work_dir, data,
                        headers=("#",)):
    """Limit CNVkit input to calls with support from another caller.

    get_coords is a function that return chrom, start, end from a line of the
    input_file, allowing handling of multiple input file types.
    """
    support_files = [(c, tz.get_in([c, "vrn_file"], calls_by_name))
                     for c in convert.SUBSET_BY_SUPPORT["cnvkit"]]
    support_files = [(c, f) for (c, f) in support_files if f and vcfutils.vcf_has_variants(f)]
    if len(support_files) == 0:
        return input_file
    else:
        out_file = os.path.join(work_dir, "%s-havesupport%s" %
                                utils.splitext_plus(os.path.basename(input_file)))
        if not utils.file_uptodate(out_file, input_file):
            input_bed = _input_to_bed(input_file, work_dir, get_coords, headers)
            pass_coords = set([])
            with file_transaction(data, out_file) as tx_out_file:
                support_beds = " ".join([_sv_vcf_to_bed(f, c, out_file) for c, f in support_files])
                tmp_cmp_bed = "%s-intersectwith.bed" % utils.splitext_plus(tx_out_file)[0]
                cmd = "bedtools intersect -wa -f 0.5 -r -a {input_bed} -b {support_beds} > {tmp_cmp_bed}"
                do.run(cmd.format(**locals()), "Intersect CNVs with support files")
                for r in pybedtools.BedTool(tmp_cmp_bed):
                    pass_coords.add((str(r.chrom), str(r.start), str(r.stop)))
                with open(input_file) as in_handle:
                    with open(tx_out_file, "w") as out_handle:
                        for line in in_handle:
                            passes = True
                            if not line.startswith(headers):
                                passes = get_coords(line) in pass_coords
                            if passes:
                                out_handle.write(line)
        return out_file

def _sv_vcf_to_bed(orig_vcf, caller, base_out):
    out_file = "%s-inputcmp-%s.bed" % (utils.splitext_plus(base_out)[0], caller)
    if not utils.file_exists(out_file):
        with pysam.VariantFile(orig_vcf) as bcf_in:
            with open(out_file, "w") as out_handle:
                for rec in bcf_in:
                    if len(rec.filter.keys()) == 0 or rec.filter.keys()[0] in ["PASS", "."]:
                        out_handle.write("%s\t%s\t%s\n" % (rec.chrom, rec.start, rec.info.get("END", rec.start)))
    return out_file

def _input_to_bed(theta_input, work_dir, get_coords, headers):
    """Convert input file to a BED file for comparisons
    """
    theta_bed = os.path.join(work_dir, "%s.bed" % os.path.splitext(os.path.basename(theta_input))[0])
    with open(theta_input) as in_handle:
        with open(theta_bed, "w") as out_handle:
            for line in in_handle:
                if not line.startswith(headers):
                    chrom, start, end = get_coords(line)
                    out_handle.write("\t".join([chrom, start, end]) + "\n")
    return theta_bed

def _run_theta(cnv_info, data, work_dir, run_n3=True):
    """Run theta, calculating subpopulations and normal contamination.
    """
    out = {"caller": "theta"}
    max_normal = "0.9"
    opts = ["-m", max_normal]
    n2_result = _safe_run_theta(cnv_info["theta_input"], os.path.join(work_dir, "n2"), ".n2.results",
                                ["-n", "2"] + opts, data)
    if n2_result:
        out["estimate"] = n2_result
        if run_n3:
            n2_bounds = "%s.withBounds" % os.path.splitext(n2_result)[0]
            n3_result = _safe_run_theta(n2_bounds, os.path.join(work_dir, "n3"), ".n3.results",
                                        ["-n", "3", "--RESULTS", n2_result] + opts,
                                        data)
            if n3_result:
                best_result = _select_model(n2_bounds, n2_result, n3_result,
                                            os.path.join(work_dir, "n3"), data)
                out["estimate"] = best_result
                out["cnvs"] = _merge_theta_calls(n2_bounds, best_result, cnv_info["vrn_file"], data)
    return out

def _update_with_calls(result_file, cnv_file):
    """Update bounds with calls from CNVkit, inferred copy numbers and p-values from THetA.
    """
    results = {}
    with open(result_file) as in_handle:
        in_handle.readline()  # header
        _, _, cs, ps = in_handle.readline().strip().split()
        for i, (c, p) in enumerate(zip(cs.split(":"), ps.split(","))):
            results[i] = (c, p)
    cnvs = {}
    with open(cnv_file) as in_handle:
        for line in in_handle:
            chrom, start, end, _, count = line.rstrip().split()[:5]
            cnvs[(chrom, start, end)] = count
    def update(i, line):
        parts = line.rstrip().split("\t")
        chrom, start, end = parts[1:4]
        parts += cnvs.get((chrom, start, end), ".")
        parts += list(results[i])
        return "\t".join(parts) + "\n"
    return update

def _merge_theta_calls(bounds_file, result_file, cnv_file, data):
    """Create a final output file with merged CNVkit and THetA copy and population estimates.
    """
    out_file = "%s-merged.txt" % (result_file.replace(".BEST.results", ""))
    if not utils.file_uptodate(out_file, result_file):
        with file_transaction(data, out_file) as tx_out_file:
            updater = _update_with_calls(result_file, cnv_file)
            with open(bounds_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    i = 0
                    for line in in_handle:
                        if line.startswith("#"):
                            parts = line.rstrip().split("\t")
                            parts += ["cnv", "pop_cnvs", "pop_pvals"]
                            out_handle.write("\t".join(parts) + "\n")
                        else:
                            out_handle.write(updater(i, line))
                            i += 1
    return out_file

def _select_model(n2_bounds, n2_result, n3_result, out_dir, data):
    """Run final model selection from n=2 and n=3 options.
    """
    n2_out_file = n2_result.replace(".n2.results", ".BEST.results")
    n3_out_file = n3_result.replace(".n3.results", ".BEST.results")
    if not utils.file_exists(n2_out_file) and not utils.file_exists(n3_out_file):
        cmd = _get_cmd("ModelSelection.py") + [n2_bounds, n2_result, n3_result]
        do.run(cmd, "Select best THetA model")
    if utils.file_exists(n2_out_file):
        return n2_out_file
    else:
        assert utils.file_exists(n3_out_file)
        return n3_out_file

def _safe_run_theta(input_file, out_dir, output_ext, args, data):
    """Run THetA, catching and continuing on any errors.
    """
    out_file = os.path.join(out_dir, _split_theta_ext(input_file) + output_ext)
    skip_file = out_file + ".skipped"
    if utils.file_exists(skip_file):
        return None
    if not utils.file_exists(out_file):
        with file_transaction(data, out_dir) as tx_out_dir:
            utils.safe_makedir(tx_out_dir)
            cmd = _get_cmd("RunTHetA.py") + args + \
                  [input_file, "--NUM_PROCESSES", dd.get_cores(data),
                   "--FORCE", "-d", tx_out_dir]
            try:
                do.run(cmd, "Run THetA to calculate purity", log_error=False)
            except subprocess.CalledProcessError, msg:
                if ("Number of intervals must be greater than 1" in str(msg) or
                      "This sample isn't a good candidate for THetA analysis" in str(msg)):
                    with open(os.path.join(tx_out_dir, os.path.basename(skip_file)), "w") as out_handle:
                        out_handle.write("Expected TheTA failure, skipping")
                    return None
                else:
                    raise
    return out_file

def _split_theta_ext(fname):
    base = os.path.splitext(os.path.basename(fname))[0]
    if base.endswith((".n2", ".n3")):
        base = os.path.splitext(base)[0]
    return base

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "heterogeneity",
                                           dd.get_sample_name(data), "theta"))

def _get_cmd(cmd):
    """Retrieve required commands for running THetA with our local bcbio python.
    """
    check_cmd = "RunTHetA.py"
    try:
        local_cmd = subprocess.check_output(["which", check_cmd]).strip()
    except subprocess.CalledProcessError:
        return None
    return [sys.executable, "%s/%s" % (os.path.dirname(os.path.realpath(local_cmd)), cmd)]
