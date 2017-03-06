"""Copy number detection with CNVkit with specific support for targeted sequencing.

http://cnvkit.readthedocs.org
"""
import copy
import math
import operator
import os
import sys
import tempfile
import subprocess

import pybedtools
import numpy as np
import toolz as tz

from bcbio import utils
from bcbio.bam import ref
from bcbio.distributed.multi import run_multicore, zeromq_aware_logging
from bcbio.distributed.transaction import file_transaction
from bcbio.heterogeneity import chromhacks
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.variation import bedutils, effects, ploidy, population, vcfutils
from bcbio.structural import annotate, shared, plot

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

def _associate_cnvkit_out(ckouts, items, is_somatic=False):
    """Associate cnvkit output with individual items.
    """
    assert len(ckouts) == len(items)
    out = []
    for ckout, data in zip(ckouts, items):
        ckout = copy.deepcopy(ckout)
        ckout["variantcaller"] = "cnvkit"
        if utils.file_exists(ckout["cns"]) and _cna_has_values(ckout["cns"]):
            ckout = _add_seg_to_output(ckout, data)
            ckout = _add_gainloss_to_output(ckout, data)
            ckout = _add_segmetrics_to_output(ckout, data)
            ckout = _add_variantcalls_to_output(ckout, data, is_somatic)
            # ckout = _add_coverage_bedgraph_to_output(ckout, data)
            ckout = _add_cnr_bedgraph_and_bed_to_output(ckout, data)
            if "svplots" in dd.get_tools_on(data):
                ckout = _add_plots_to_output(ckout, data)
            if "sv" not in data:
                data["sv"] = []
            data["sv"].append(ckout)
        out.append(data)
    return out

def _run_cnvkit_single(data, background=None):
    """Process a single input file with BAM or uniform background.
    """
    if not background:
        background = []
    ckouts = _run_cnvkit_shared([data], background)
    if not ckouts:
        return [data]
    else:
        assert len(ckouts) == 1
        return _associate_cnvkit_out(ckouts, [data])

def _run_cnvkit_cancer(items, background):
    """Run CNVkit on a tumor/normal pair.
    """
    paired = vcfutils.get_paired_bams([x["align_bam"] for x in items], items)
    normal_data = [x for x in items if dd.get_sample_name(x) != paired.tumor_name]
    tumor_ready, normal_ready = _match_batches(paired.tumor_data, normal_data[0] if normal_data else None)
    ckouts = _run_cnvkit_shared([tumor_ready], [normal_ready] if normal_ready else [])
    if not ckouts:
        return items
    assert len(ckouts) == 1
    tumor_data = _associate_cnvkit_out(ckouts, [paired.tumor_data], is_somatic=True)
    return tumor_data + normal_data

def _match_batches(tumor, normal):
    """Fix batch names for shared tumor/normals to ensure matching
     """
    def _get_batch(x):
        b = dd.get_batch(x)
        return [b] if not isinstance(b, (list, tuple)) else b
    if normal:
        tumor = copy.deepcopy(tumor)
        normal = copy.deepcopy(normal)
        cur_batch = list(set(_get_batch(tumor)) & set(_get_batch(normal)))
        assert len(cur_batch) == 1, "No batch overlap: %s and %s" % (_get_batch(tumor), _get_batch(normal))
        cur_batch = cur_batch[0]
        tumor["metadata"]["batch"] = cur_batch
        normal["metadata"]["batch"] = cur_batch
    return tumor, normal

def _run_cnvkit_population(items, background):
    """Run CNVkit on a population of samples.

    Tries to calculate background based on case/controls, otherwise
    uses samples from the same batch as background.
    """
    if background and len(background) > 0:
        inputs = items
    else:
        inputs, background = shared.find_case_control(items)

    # if we have case/control organized background or a single sample
    if len(inputs) == 1 or len(background) > 0:
        ckouts = _run_cnvkit_shared(inputs, background)
        return _associate_cnvkit_out(ckouts, inputs) + background
    # otherwise run each sample with the others in the batch as background
    else:
        out = []
        for cur_input in items:
            background = [d for d in items if dd.get_sample_name(d) != dd.get_sample_name(cur_input)]
            ckouts = _run_cnvkit_shared([cur_input], background)
            out.extend(_associate_cnvkit_out(ckouts, [cur_input]))
        return out

def _get_cmd(script_name="cnvkit.py"):
    return os.path.join(os.path.dirname(os.path.realpath(sys.executable)), script_name)

def _prep_cmd(cmd, tx_out_file):
    """Wrap CNVkit commands ensuring we use local temporary directories.
    """
    cmd = " ".join(cmd) if isinstance(cmd, (list, tuple)) else cmd
    return "export TMPDIR=%s && %s" % (os.path.dirname(tx_out_file), cmd)

def _bam_to_outbase(bam_file, work_dir):
    """Convert an input BAM file into CNVkit expected output.
    """
    out_base = os.path.splitext(os.path.basename(bam_file))[0].split(".")[0]
    return os.path.join(work_dir, out_base)

def _run_cnvkit_shared(inputs, backgrounds):
    """Shared functionality to run CNVkit, parallelizing over multiple BAM files.
    """
    work_dir = _sv_workdir(inputs[0])
    raw_work_dir = utils.safe_makedir(os.path.join(work_dir, "raw"))
    background_name = dd.get_sample_name(backgrounds[0]) if backgrounds else "flat"
    background_cnn = os.path.join(raw_work_dir, "%s_background.cnn" % (background_name))
    ckouts = []
    for cur_input in inputs:
        cur_raw_work_dir = utils.safe_makedir(os.path.join(_sv_workdir(cur_input), "raw"))
        out_base = _bam_to_outbase(dd.get_align_bam(cur_input), cur_raw_work_dir)
        ckouts.append({"cnr": "%s.cnr" % out_base,
                       "cns": "%s.cns" % out_base,
                       "back_cnn": background_cnn})
    if not utils.file_exists(ckouts[0]["cns"]):
        cov_interval = dd.get_coverage_interval(inputs[0])
        raw_target_bed, access_bed = _get_target_access_files(cov_interval, inputs[0], work_dir)
        # bail out if we ended up with no regions
        if not utils.file_exists(raw_target_bed):
            return {}
        raw_target_bed = annotate.add_genes(raw_target_bed, inputs[0])
        parallel = {"type": "local", "cores": dd.get_cores(inputs[0]), "progs": ["cnvkit"]}
        target_bed, antitarget_bed = _cnvkit_targets(raw_target_bed, access_bed, cov_interval,
                                                     raw_work_dir, inputs[0])
        samples_to_run = zip(["background"] * len(backgrounds), backgrounds) + \
                         zip(["evaluate"] * len(inputs), inputs)
        raw_coverage_cnns = [_cnvkit_coverage(cdata, bed, itype) for itype, cdata in samples_to_run
                             for bed in [target_bed, antitarget_bed]]
        coverage_cnns = reduce(operator.add,
                               [_cnvkit_metrics(cnns, target_bed, antitarget_bed, cov_interval, inputs + backgrounds)
                                for cnns in tz.groupby("bam", raw_coverage_cnns).values()])
        background_cnn = _cnvkit_background(_select_background_cnns(coverage_cnns),
                                            background_cnn, target_bed, antitarget_bed, inputs[0])
        fixed_cnrs = run_multicore(_cnvkit_fix,
                                   [(cnns, background_cnn, inputs + backgrounds) for cnns in
                                    tz.groupby("bam", [x for x in coverage_cnns
                                                       if x["itype"] == "evaluate"]).values()],
                                      inputs[0]["config"], parallel)
        [_cnvkit_segment(cnr, cov_interval, data) for cnr, data in fixed_cnrs]
    return ckouts

def _cna_has_values(fname):
    with open(fname) as in_handle:
        for i, line in enumerate(in_handle):
            if i > 0:
                return True
    return False

def _cnvkit_segment(cnr_file, cov_interval, data):
    """Perform segmentation and copy number calling on normalized inputs
    """
    out_file = "%s.cns" % os.path.splitext(cnr_file)[0]
    if not utils.file_uptodate(out_file, cnr_file):
        with file_transaction(data, out_file) as tx_out_file:
            if not _cna_has_values(cnr_file):
                with open(tx_out_file, "w") as out_handle:
                    out_handle.write("chromosome\tstart\tend\tgene\tlog2\tprobes\tCN1\tCN2\tbaf\tweight\n")
            else:
                cmd = [_get_cmd(), "segment", "-p", str(dd.get_cores(data)),
                       "-o", tx_out_file, cnr_file]
                small_vrn_files = _compatible_small_variants(data)
                if len(small_vrn_files) > 0 and _cna_has_values(cnr_file) and cov_interval != "genome":
                    cmd += ["-v", small_vrn_files[0]]
                if cov_interval == "genome":
                    cmd += ["--threshold", "0.00001"]
                # preferentially use conda installed Rscript
                export_cmd = ("unset R_HOME && export PATH=%s:$PATH && export TMPDIR=%s && "
                              % (os.path.dirname(utils.Rscript_cmd()), os.path.dirname(tx_out_file)))
                do.run(export_cmd + " ".join(cmd), "CNVkit segment")
    return out_file

def _cnvkit_metrics(cnns, target_bed, antitarget_bed, cov_interval, items):
    """Estimate noise of a sample using a flat background.

    Only used for panel/targeted data due to memory issues with whole genome
    samples.
    """
    if cov_interval == "genome":
        return cnns
    target_cnn = [x["file"] for x in cnns if x["cnntype"] == "target"][0]
    background_file = "%s-flatbackground.cnn" % utils.splitext_plus(target_cnn)[0]
    background_file = _cnvkit_background([], background_file, target_bed, antitarget_bed, items[0])
    cnr_file, data = _cnvkit_fix_base(cnns, background_file, items, "-flatbackground")
    cns_file = _cnvkit_segment(cnr_file, cov_interval, data)
    metrics_file = "%s-metrics.txt" % utils.splitext_plus(target_cnn)[0]
    if not utils.file_exists(metrics_file):
        with file_transaction(data, metrics_file) as tx_metrics_file:
            cmd = [_get_cmd(), "metrics", "-o", tx_metrics_file, "-s", cns_file, "--", cnr_file]
            do.run(_prep_cmd(cmd, tx_metrics_file), "CNVkit metrics")
    metrics = _read_metrics_file(metrics_file)
    out = []
    for cnn in cnns:
        cnn["metrics"] = metrics
        out.append(cnn)
    return out

def _read_metrics_file(in_file):
    with open(in_file) as in_handle:
        header = in_handle.next().strip().split("\t")[1:]
        vals = map(float, in_handle.next().strip().split("\t")[1:])
    return dict(zip(header, vals))

@utils.map_wrap
@zeromq_aware_logging
def _cnvkit_fix(cnns, background_cnn, items):
    """Normalize samples, correcting sources of bias.
    """
    return [_cnvkit_fix_base(cnns, background_cnn, items)]

def _cnvkit_fix_base(cnns, background_cnn, items, ext=""):
    assert len(cnns) == 2, "Expected target and antitarget CNNs: %s" % cnns
    target_cnn = [x["file"] for x in cnns if x["cnntype"] == "target"][0]
    antitarget_cnn = [x["file"] for x in cnns if x["cnntype"] == "antitarget"][0]
    data = [x for x in items if dd.get_sample_name(x) == cnns[0]["sample"]][0]
    out_file = "%s%s.cnr" % (os.path.commonprefix([target_cnn, antitarget_cnn]).replace(".", ""), ext)
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [_get_cmd(), "fix", "-o", tx_out_file, target_cnn, antitarget_cnn, background_cnn]
            do.run(_prep_cmd(cmd, tx_out_file), "CNVkit fix")
    return out_file, data

def _select_background_cnns(cnns):
    """Select cnns to use for background calculations.

    Uses background samples in cohort, and will remove CNNs with high
    on target variability. Uses (number of segments * biweight midvariance) as metric
    for variability with higher numbers being more unreliable.
    """
    min_for_variability_analysis = 20
    pct_keep = 0.10
    b_cnns = [x for x in cnns if x["itype"] == "background" and x.get("metrics")]
    assert len(b_cnns) % 2 == 0, "Expect even set of target/antitarget cnns for background"
    if len(b_cnns) >= min_for_variability_analysis:
        b_cnns_w_metrics = []
        for b_cnn in b_cnns:
            unreliability = b_cnn["metrics"]["segments"] * b_cnn["metrics"]["bivar"]
            b_cnns_w_metrics.append((unreliability, b_cnn))
        b_cnns_w_metrics.sort()
        to_keep = int(math.ceil(pct_keep * len(b_cnns) / 2.0) * 2)
        b_cnns = [x[1] for x in b_cnns_w_metrics][:to_keep]
        assert len(b_cnns) % 2 == 0, "Expect even set of target/antitarget cnns for background"
    return [x["file"] for x in b_cnns]

def _cnvkit_background(background_cnns, out_file, target_bed, antitarget_bed, data):
    """Calculate background reference, handling flat case with no normal sample.
    """
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [_get_cmd(), "reference", "-f", dd.get_ref_file(data), "-o", tx_out_file]
            if len(background_cnns) == 0:
                cmd += ["-t", target_bed, "-a", antitarget_bed]
            else:
                cmd += background_cnns
            do.run(_prep_cmd(cmd, tx_out_file), "CNVkit background")
    return out_file

def _cnvkit_coverage(data, bed_file, input_type):
    """Calculate coverage in a BED file for CNVkit.
    """
    bam_file = dd.get_align_bam(data)
    work_dir = utils.safe_makedir(os.path.join(_sv_workdir(data), "raw"))
    exts = {".target.bed": ("target", "targetcoverage.cnn"),
            ".antitarget.bed": ("antitarget", "antitargetcoverage.cnn")}
    cnntype = None
    for orig, (cur_cnntype, ext) in exts.items():
        if bed_file.endswith(orig):
            cnntype = cur_cnntype
            break
    if cnntype is None:
        assert bed_file.endswith(".bed"), "Unexpected BED file extension for coverage %s" % bed_file
        cnntype = ""
    batch = dd.get_batch(data) or dd.get_sample_name(data)
    base = _bam_to_outbase(bam_file, work_dir)
    out_file = "%s-%s.%s" % (base, batch, ext)
    out_file_old = "%s.%s" % (base, ext)
    # back compatible with previous runs to avoid re-calculating
    if utils.file_exists(out_file_old):
        out_file = out_file_old
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [_get_cmd(), "coverage", "-p", str(dd.get_cores(data)), bam_file, bed_file, "-o", tx_out_file]
            do.run(_prep_cmd(cmd, tx_out_file), "CNVkit coverage")
    return {"itype": input_type, "file": out_file, "bam": bam_file, "cnntype": cnntype,
            "sample": dd.get_sample_name(data)}

def _cnvkit_targets(raw_target_bed, access_bed, cov_interval, work_dir, data):
    """Create target and antitarget regions from target and access files.
    """
    batch = dd.get_batch(data) or dd.get_sample_name(data)
    basename = os.path.splitext(os.path.basename(raw_target_bed))[0]
    target_bed = os.path.join(work_dir, "%s-%s.target.bed" % (basename, batch))
    # back compatible with previous runs to avoid re-calculating
    target_bed_old = os.path.join(work_dir, "%s.target.bed" % basename)
    if utils.file_exists(target_bed_old):
        target_bed = target_bed_old
    if not utils.file_exists(target_bed):
        with file_transaction(data, target_bed) as tx_out_file:
            cmd = [_get_cmd(), "target", raw_target_bed, "--split", "-o", tx_out_file]
            bin_estimates = _cnvkit_coverage_bin_estimate(raw_target_bed, access_bed, cov_interval, work_dir, data)
            if bin_estimates.get("target"):
                cmd += ["--avg-size", str(bin_estimates["target"])]
            do.run(_prep_cmd(cmd, tx_out_file), "CNVkit target")
    antitarget_bed = os.path.join(work_dir, "%s-%s.antitarget.bed" % (basename, batch))
    antitarget_bed_old = os.path.join(work_dir, "%s.antitarget.bed" % basename)
    # back compatible with previous runs to avoid re-calculating
    if os.path.exists(antitarget_bed_old):
        antitarget_bed = antitarget_bed_old
    if not os.path.exists(antitarget_bed):
        with file_transaction(data, antitarget_bed) as tx_out_file:
            cmd = [_get_cmd(), "antitarget", "-g", access_bed, target_bed, "-o", tx_out_file]
            bin_estimates = _cnvkit_coverage_bin_estimate(raw_target_bed, access_bed, cov_interval, work_dir, data)
            if bin_estimates.get("antitarget"):
                cmd += ["--avg-size", str(bin_estimates["antitarget"])]
            do.run(_prep_cmd(cmd, tx_out_file), "CNVkit antitarget")
    return target_bed, antitarget_bed

def _cnvkit_coverage_bin_estimate(raw_target_bed, access_bed, cov_interval, work_dir, data):
    """Estimate good coverage bin sizes for target regions based on coverage.
    """
    batch = dd.get_batch(data) or dd.get_sample_name(data)
    out_file = os.path.join(work_dir, "%s-%s-bin_estimate.txt" % (
        os.path.splitext(os.path.basename(raw_target_bed))[0], batch))
    method_map = {"genome": "wgs", "regional": "hybrid", "amplicon": "amplicon"}
    if not os.path.exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [_get_cmd("coverage_bin_size.py"), dd.get_align_bam(data),
                   "-m", method_map[cov_interval], "-t", raw_target_bed,
                   "-g", access_bed]
            cmd = " ".join(cmd) + " > " + tx_out_file
            try:
                do.run(_prep_cmd(cmd, tx_out_file), "CNVkit coverage bin estimation", log_error=False)
            except subprocess.CalledProcessError:
                logger.info("Bin size estimate failed, using default values")
                with open(tx_out_file, "w") as out_handle:
                    out_handle.write("Bin size estimate failed, using default values")
    avg_bin_sizes = {}
    estimate_map = {"On-target": "target", "Off-target": "antitarget",
                    "Genome": "target", "Targets (sampling)": "target"}
    range_map = {("genome", "target"): (500, 1000),
                 ("regional", "target"): (50, 267), ("regional", "antitarget"): (20000, 200000),
                 ("amplicon", "target"): (50, 267)}
    with open(out_file) as in_handle:
        for line in in_handle:
            if line.startswith(tuple(estimate_map.keys())):
                name, depth, bin_size = line.strip().split("\t")
                name = estimate_map[name.replace(":", "").strip()]
                try:
                    bin_size = int(bin_size)
                except ValueError:
                    bin_size = None
                if bin_size and bin_size > 0:
                    cur_min, cur_max = range_map[(cov_interval, name)]
                    avg_bin_sizes[name] = max(min(bin_size, cur_max), cur_min)
    return avg_bin_sizes

def _get_target_access_files(cov_interval, data, work_dir):
    """Retrieve target and access files based on the type of data to process.

    pick targets, anti-targets and access files based on analysis type
    http://cnvkit.readthedocs.org/en/latest/nonhybrid.html
    """
    base_regions = shared.get_base_cnv_regions(data, work_dir)
    target_bed = bedutils.sort_merge(base_regions, data, out_dir=work_dir)
    if cov_interval == "amplicon":
        return target_bed, target_bed
    elif cov_interval == "genome":
        return target_bed, target_bed
    else:
        access_file = _create_access_file(dd.get_ref_file(data), _sv_workdir(data), data)
        return target_bed, access_file

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

def _add_cnr_bedgraph_and_bed_to_output(out, data):
    cnr_file = out["cnr"]
    bedgraph_file = cnr_file + ".bedgraph"
    if not utils.file_exists(bedgraph_file):
        with file_transaction(data, bedgraph_file) as tx_out_file:
            cmd = "sed 1d {cnr_file} | cut -f1,2,3,5 > {tx_out_file}"
            do.run(cmd.format(**locals()), "Converting cnr to bedgraph format")
    out["cnr_bedgraph"] = bedgraph_file

    bed_file = cnr_file + ".bed"
    if not utils.file_exists(bed_file):
        with file_transaction(data, bed_file) as tx_out_file:
            cmd = "sed 1d {cnr_file} | cut -f1,2,3,4,5 > {tx_out_file}"
            do.run(cmd.format(**locals()), "Converting cnr to bed format")
    out["cnr_bed"] = bed_file
    return out

def _compatible_small_variants(data):
    """Retrieve small variant (SNP, indel) VCFs compatible with CNVkit.
    """
    supported = set(["vardict", "freebayes", "gatk-haplotype", "mutect2", "vardict"])
    out = []
    for v in data.get("variants", []):
        vrn_file = v.get("vrn_file")
        if vrn_file and v.get("variantcaller") in supported:
            base, ext = utils.splitext_plus(os.path.basename(vrn_file))
            if vcfutils.get_paired_phenotype(data):
                out.append(vrn_file)
            else:
                sample_vrn_file = os.path.join(dd.get_work_dir(data), v["variantcaller"],
                                               "%s-%s%s" % (base, dd.get_sample_name(data), ext))
                sample_vrn_file = vcfutils.select_sample(vrn_file, dd.get_sample_name(data), sample_vrn_file,
                                                         data["config"])
                out.append(sample_vrn_file)
    return out

def _add_variantcalls_to_output(out, data, is_somatic=False):
    """Call ploidy and convert into VCF and BED representations.
    """
    call_file = "%s-call%s" % os.path.splitext(out["cns"])
    gender = population.get_gender(data)
    if not utils.file_exists(call_file):
        with file_transaction(data, call_file) as tx_call_file:
            cmd = [os.path.join(os.path.dirname(sys.executable), "cnvkit.py"), "call",
                   "--ploidy", str(ploidy.get_ploidy([data])),
                   "-o", tx_call_file, out["cns"]]
            small_vrn_files = _compatible_small_variants(data)
            if len(small_vrn_files) > 0 and _cna_has_values(out["cns"]):
                cmd += ["-v", small_vrn_files[0]]
                if not is_somatic:
                    cmd += ["-m", "clonal"]
            if gender and gender.lower() != "unknown":
                cmd += ["--gender", gender]
                if gender.lower() == "male":
                    cmd += ["--male-reference"]
            do.run(cmd, "CNVkit call ploidy")
    calls = {}
    for outformat in ["bed", "vcf"]:
        out_file = "%s.%s" % (os.path.splitext(call_file)[0], outformat)
        calls[outformat] = out_file
        if not os.path.exists(out_file):
            with file_transaction(data, out_file) as tx_out_file:
                cmd = [os.path.join(os.path.dirname(sys.executable), "cnvkit.py"), "export",
                       outformat, "--sample-id", dd.get_sample_name(data),
                       "--ploidy", str(ploidy.get_ploidy([data])),
                       "-o", tx_out_file, call_file]
                if gender and gender.lower() == "male":
                    cmd += ["--male-reference"]
                do.run(cmd, "CNVkit export %s" % outformat)
    out["call_file"] = call_file
    out["vrn_bed"] = annotate.add_genes(calls["bed"], data)
    effects_vcf, _ = effects.add_to_vcf(calls["vcf"], data, "snpeff")
    out["vrn_file"] = effects_vcf or calls["vcf"]
    return out

def _add_segmetrics_to_output(out, data):
    """Add metrics for measuring reliability of CNV estimates.
    """
    out_file = "%s-segmetrics.txt" % os.path.splitext(out["cns"])[0]
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [os.path.join(os.path.dirname(sys.executable), "cnvkit.py"), "segmetrics",
                   "--ci", "--pi",
                   "-s", out["cns"], "-o", tx_out_file, out["cnr"]]
            # Use less fine grained bootstrapping intervals for whole genome runs
            if dd.get_coverage_interval(data) == "genome":
                cmd += ["--alpha", "0.1", "--bootstrap", "50"]
            else:
                cmd += ["--alpha", "0.01", "--bootstrap", "500"]
            do.run(cmd, "CNVkit segmetrics")
    out["segmetrics"] = out_file
    return out

def _add_gainloss_to_output(out, data):
    """Add gainloss based on genes, helpful for identifying changes in smaller genes.
    """
    out_file = "%s-gainloss.txt" % os.path.splitext(out["cns"])[0]
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [os.path.join(os.path.dirname(sys.executable), "cnvkit.py"), "gainloss",
                   "-s", out["cns"], "-o", tx_out_file, out["cnr"]]
            do.run(cmd, "CNVkit gainloss")
    out["gainloss"] = out_file
    return out

def _add_coverage_bedgraph_to_output(out, data):
    """Add BedGraph representation of coverage to the output
    """
    out_file = "%s.coverage.bedgraph" % os.path.splitext(out["cns"])[0]
    if utils.file_exists(out_file):
        out["bedgraph"] = out_file
        return out
    bam_file = dd.get_align_bam(data)
    bedtools = config_utils.get_program("bedtools", data["config"])
    samtools = config_utils.get_program("samtools", data["config"])
    cns_file = out["cns"]
    bed_file = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
    with file_transaction(data, out_file) as tx_out_file:
        cmd = ("sed 1d {cns_file} | cut -f1,2,3 > {bed_file}; "
               "{samtools} view -b -L {bed_file} {bam_file} | "
               "{bedtools} genomecov -bg -ibam - -g {bed_file} >"
               "{tx_out_file}").format(**locals())
        do.run(cmd, "CNVkit bedGraph conversion")
        os.remove(bed_file)
    out["bedgraph"] = out_file
    return out

def _add_plots_to_output(out, data):
    """Add CNVkit plots summarizing called copy number values.
    """
    out["plot"] = {}
    diagram_plot = _add_diagram_plot(out, data)
    if diagram_plot:
        out["plot"]["diagram"] = diagram_plot
    scatter = _add_scatter_plot(out, data)
    if scatter:
        out["plot"]["scatter"] = scatter
    scatter_global = _add_global_scatter_plot(out, data)
    if scatter_global:
        out["plot"]["scatter_global"] = scatter_global
    return out

def _get_larger_chroms(ref_file):
    """Retrieve larger chromosomes, avoiding the smaller ones for plotting.
    """
    from scipy.cluster.vq import kmeans, vq
    all_sizes = []
    for c in ref.file_contigs(ref_file):
        all_sizes.append(float(c.size))
    all_sizes.sort()
    # separate out smaller chromosomes and haplotypes with kmeans
    centroids, _ = kmeans(np.array(all_sizes), 2)
    idx, _ = vq(np.array(all_sizes), centroids)
    little_sizes = tz.first(tz.partitionby(lambda xs: xs[0], zip(idx, all_sizes)))
    little_sizes = [x[1] for x in little_sizes]
    # create one more cluster with the smaller, removing the haplotypes
    centroids2, _ = kmeans(np.array(little_sizes), 2)
    idx2, _ = vq(np.array(little_sizes), centroids2)
    little_sizes2 = tz.first(tz.partitionby(lambda xs: xs[0], zip(idx2, little_sizes)))
    little_sizes2 = [x[1] for x in little_sizes2]
    # get any chromosomes not in haplotype/random bin
    thresh = max(little_sizes2)
    larger_chroms = []
    for c in ref.file_contigs(ref_file):
        if c.size > thresh:
            larger_chroms.append(c.name)
    return larger_chroms

def _remove_haplotype_chroms(in_file, data):
    """Remove shorter haplotype chromosomes from cns/cnr files for plotting.
    """
    larger_chroms = set(_get_larger_chroms(dd.get_ref_file(data)))
    out_file = "%s-chromfilter%s" % utils.splitext_plus(in_file)
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(in_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    for line in in_handle:
                        if line.startswith("chromosome") or line.split()[0] in larger_chroms:
                            out_handle.write(line)
    return out_file

def _add_global_scatter_plot(out, data):
    out_file = "%s-scatter_global.pdf" % os.path.splitext(out["cnr"])[0]
    if utils.file_exists(out_file):
        return out_file
    cnr = _remove_haplotype_chroms(out["cnr"], data)
    cns = _remove_haplotype_chroms(out["cns"], data)
    with file_transaction(data, out_file) as tx_out_file:
        cmd = [_get_cmd(), "scatter", "-s", cns, "-o", tx_out_file, cnr]
        do.run(_prep_cmd(cmd, tx_out_file), "CNVkit global scatter plot")
    return out_file

def _add_scatter_plot(out, data):
    out_file = "%s-scatter.pdf" % os.path.splitext(out["cnr"])[0]
    priority_bed = dd.get_svprioritize(data)
    if not priority_bed:
        return None
    priority_bed = plot._prioritize_plot_regions(pybedtools.BedTool(priority_bed), data, os.path.dirname(out_file))
    if utils.file_exists(out_file):
        return out_file
    cnr = _remove_haplotype_chroms(out["cnr"], data)
    cns = _remove_haplotype_chroms(out["cns"], data)
    with file_transaction(data, out_file) as tx_out_file:
        cmd = [_get_cmd(), "scatter", "-s", cns, "-o", tx_out_file, "-l",
               priority_bed, cnr]
        do.run(_prep_cmd(cmd, tx_out_file), "CNVkit scatter plot")
    return out_file

def _cnx_is_empty(in_file):
    """Check if cnr or cns files are empty (only have a header)
    """
    with open(in_file) as in_handle:
        for i, line in enumerate(in_handle):
            if i > 0:
                return False
    return True

def _add_diagram_plot(out, data):
    out_file = "%s-diagram.pdf" % os.path.splitext(out["cnr"])[0]
    cnr = _remove_haplotype_chroms(out["cnr"], data)
    cns = _remove_haplotype_chroms(out["cns"], data)
    if _cnx_is_empty(cnr) or _cnx_is_empty(cns):
        return None
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [_get_cmd(), "diagram", "-s", cns,
                   "-o", tx_out_file, cnr]
            gender = population.get_gender(data)
            if gender and gender.lower() == "male":
                cmd += ["--male-reference"]
            do.run(_prep_cmd(cmd, tx_out_file), "CNVkit diagram plot")
    return out_file

def _create_access_file(ref_file, out_dir, data):
    """Create genome access file for CNVlib to define available genomic regions.

    XXX Can move to installation/upgrade process if too slow here.
    """
    out_file = os.path.join(out_dir, "%s-access.bed" % os.path.splitext(os.path.basename(ref_file))[0])
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [_get_cmd(), "access",
                   ref_file, "-s", "10000", "-o", tx_out_file]
            do.run(_prep_cmd(cmd, tx_out_file), "Create CNVkit access file")
    return out_file

# ## Theta support

def export_theta(ckout, data):
    """Provide updated set of data with export information for TheTA2 input.
    """
    cns_file = chromhacks.bed_to_standardonly(ckout["cns"], data, headers="chromosome")
    cnr_file = chromhacks.bed_to_standardonly(ckout["cnr"], data, headers="chromosome")
    out_file = "%s-theta.input" % utils.splitext_plus(cns_file)[0]
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [_get_cmd(), "export", "theta", cns_file, cnr_file, "-o", tx_out_file]
            do.run(_prep_cmd(cmd, tx_out_file), "Export CNVkit calls as inputs for TheTA2")
    ckout["theta_input"] = out_file
    return ckout
