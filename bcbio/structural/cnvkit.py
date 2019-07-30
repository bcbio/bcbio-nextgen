"""Copy number detection with CNVkit with specific support for targeted sequencing.

http://cnvkit.readthedocs.org
"""
import collections
import copy
import glob
import math
import operator
import os
import shutil
import sys
import tempfile

import pybedtools
import numpy as np
import toolz as tz

from bcbio import utils
from bcbio.bam import ref
from bcbio.distributed.multi import run_multicore, zeromq_aware_logging
from bcbio.distributed.transaction import file_transaction
from bcbio.heterogeneity import chromhacks
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.variation import effects, ploidy, population, vcfutils
from bcbio.structural import annotate, plot, shared
from functools import reduce

def use_general_sv_bins(data):
    """Check if we should use a general binning approach for a sample.

    Checks if CNVkit is enabled and we haven't already run CNVkit.
    """
    if any([c in dd.get_svcaller(data) for c in ["cnvkit", "titancna", "purecn", "gatk-cnv"]]):
        if not _get_original_coverage(data):
            return True
    return False

def bin_approach(data):
    """Check for binning approach from configuration or normalized file.
    """
    for approach in ["cnvkit", "gatk-cnv"]:
        if approach in dd.get_svcaller(data):
            return approach
    norm_file = tz.get_in(["depth", "bins", "normalized"], data)
    if norm_file.endswith(("-crstandardized.tsv", "-crdenoised.tsv")):
        return "gatk-cnv"
    if norm_file.endswith(".cnr"):
        return "cnvkit"

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
    upload_counts = collections.defaultdict(int)
    for ckout, data in zip(ckouts, items):
        ckout = copy.deepcopy(ckout)
        ckout["variantcaller"] = "cnvkit"
        if utils.file_exists(ckout["cns"]) and _cna_has_values(ckout["cns"]):
            ckout = _add_seg_to_output(ckout, data)
            ckout = _add_gainloss_to_output(ckout, data)
            ckout = _add_segmetrics_to_output(ckout, data)
            ckout = _add_variantcalls_to_output(ckout, data, items, is_somatic)
            # ckout = _add_coverage_bedgraph_to_output(ckout, data)
            ckout = _add_cnr_bedgraph_and_bed_to_output(ckout, data)
            if "svplots" in dd.get_tools_on(data):
                ckout = _add_plots_to_output(ckout, data)
            ckout["do_upload"] = upload_counts[ckout.get("vrn_file")] == 0
        if "sv" not in data:
            data["sv"] = []
        data["sv"].append(ckout)
        if ckout.get("vrn_file"):
            upload_counts[ckout["vrn_file"]] += 1
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

def _bam_to_outbase(bam_file, work_dir, data):
    """Convert an input BAM file into CNVkit expected output.

    Handles previous non-batch cases to avoid re-calculating,
    returning both new and old values:
    """
    batch = dd.get_batch(data) or dd.get_sample_name(data)
    out_base = os.path.splitext(os.path.basename(bam_file))[0].split(".")[0]
    base = os.path.join(work_dir, out_base)
    return "%s-%s" % (base, batch), base

def _get_original_coverage(data, itype="target"):
    """Back compatible: get existing coverage files if they exist
    """
    work_dir = os.path.join(_sv_workdir(data), "raw")
    work_bam = dd.get_work_bam(data) or dd.get_align_bam(data)
    out = []
    base, _ = _bam_to_outbase(work_bam, work_dir, data)
    target_cnn = "%s.targetcoverage.cnn" % base
    anti_cnn = "%s.antitargetcoverage.cnn" % base
    if os.path.exists(target_cnn) and os.path.exists(anti_cnn):
        out.append({"bam": work_bam, "file": target_cnn, "cnntype": "target",
                    "itype": itype, "sample": dd.get_sample_name(data)})
        out.append({"bam": work_bam, "file": anti_cnn, "cnntype": "antitarget",
                    "itype": itype, "sample": dd.get_sample_name(data)})
    return out

def _run_cnvkit_shared(inputs, backgrounds):
    """Shared functionality to run CNVkit, parallelizing over multiple BAM files.

    Handles new style cases where we have pre-normalized inputs and
    old cases where we run CNVkit individually.
    """
    if tz.get_in(["depth", "bins", "normalized"], inputs[0]):
        ckouts = []
        for data in inputs:
            cnr_file = tz.get_in(["depth", "bins", "normalized"], data)
            cns_file = os.path.join(_sv_workdir(data), "%s.cns" % dd.get_sample_name(data))
            cns_file = _cnvkit_segment(cnr_file, dd.get_coverage_interval(data), data,
                                       inputs + backgrounds, cns_file)
            ckouts.append({"cnr": cnr_file, "cns": cns_file,
                           "background": tz.get_in(["depth", "bins", "background"], data)})
        return ckouts
    else:
        return _run_cnvkit_shared_orig(inputs, backgrounds)

def _get_original_targets(data):
    """Back compatible: get pre-existing target BEDs.
    """
    work_dir = os.path.join(_sv_workdir(data), "raw")
    batch = dd.get_batch(data) or dd.get_sample_name(data)
    return (glob.glob(os.path.join(work_dir, "*-%s.target.bed" % batch))[0],
            glob.glob(os.path.join(work_dir, "*-%s.antitarget.bed" % batch))[0])

def _get_general_coverage(data, itype):
    """Retrieve coverage information from new shared SV bins.
    """
    work_bam = dd.get_align_bam(data) or dd.get_work_bam(data)
    return [{"bam": work_bam, "file": tz.get_in(["depth", "bins", "target"], data),
             "cnntype": "target", "itype": itype, "sample": dd.get_sample_name(data)},
            {"bam": work_bam, "file": tz.get_in(["depth", "bins", "antitarget"], data),
             "cnntype": "antitarget", "itype": itype, "sample": dd.get_sample_name(data)}]

def _run_cnvkit_shared_orig(inputs, backgrounds):
    """Original CNVkit implementation with full normalization and segmentation.
    """
    work_dir = _sv_workdir(inputs[0])
    raw_work_dir = utils.safe_makedir(os.path.join(work_dir, "raw"))
    background_name = dd.get_sample_name(backgrounds[0]) if backgrounds else "flat"
    background_cnn = os.path.join(raw_work_dir, "%s_background.cnn" % (background_name))
    ckouts = []
    for cur_input in inputs:
        cur_raw_work_dir = utils.safe_makedir(os.path.join(_sv_workdir(cur_input), "raw"))
        out_base, out_base_old = _bam_to_outbase(dd.get_align_bam(cur_input), cur_raw_work_dir, cur_input)
        if utils.file_exists(out_base_old + ".cns"):
            out_base = out_base_old
        ckouts.append({"cnr": "%s.cnr" % out_base,
                       "cns": "%s.cns" % out_base})
    if not utils.file_exists(ckouts[0]["cns"]):
        cov_interval = dd.get_coverage_interval(inputs[0])
        samples_to_run = list(zip(["background"] * len(backgrounds), backgrounds)) + \
                         list(zip(["evaluate"] * len(inputs), inputs))
        # New style shared SV bins
        if tz.get_in(["depth", "bins", "target"], inputs[0]):
            target_bed = tz.get_in(["depth", "bins", "target"], inputs[0])
            antitarget_bed = tz.get_in(["depth", "bins", "antitarget"], inputs[0])
            raw_coverage_cnns = reduce(operator.add,
                                       [_get_general_coverage(cdata, itype) for itype, cdata in samples_to_run])
        # Back compatible with pre-existing runs
        else:
            target_bed, antitarget_bed = _get_original_targets(inputs[0])
            raw_coverage_cnns = reduce(operator.add,
                                       [_get_original_coverage(cdata, itype) for itype, cdata in samples_to_run])
        # Currently metrics not calculated due to speed and needing re-evaluation
        # We could re-enable with larger truth sets to evaluate background noise
        # But want to reimplement in a more general fashion as part of normalization
        if False:
            coverage_cnns = reduce(operator.add,
                                [_cnvkit_metrics(cnns, target_bed, antitarget_bed, cov_interval,
                                                    inputs + backgrounds)
                                    for cnns in tz.groupby("bam", raw_coverage_cnns).values()])
            background_cnn = cnvkit_background(_select_background_cnns(coverage_cnns),
                                                background_cnn, inputs, target_bed, antitarget_bed)
        else:
            coverage_cnns = raw_coverage_cnns
            background_cnn = cnvkit_background([x["file"] for x in coverage_cnns if x["itype"] == "background"],
                                                background_cnn, inputs, target_bed, antitarget_bed)
        parallel = {"type": "local", "cores": dd.get_cores(inputs[0]), "progs": ["cnvkit"]}
        fixed_cnrs = run_multicore(_cnvkit_fix,
                                   [(cnns, background_cnn, inputs, ckouts) for cnns in
                                    tz.groupby("bam", [x for x in coverage_cnns
                                                       if x["itype"] == "evaluate"]).values()],
                                   inputs[0]["config"], parallel)
        [_cnvkit_segment(cnr, cov_interval, data, inputs + backgrounds) for cnr, data in fixed_cnrs]
    return ckouts

def _cna_has_values(fname):
    with open(fname) as in_handle:
        for i, line in enumerate(in_handle):
            if i > 0:
                return True
    return False

def _cnvkit_segment(cnr_file, cov_interval, data, items, out_file=None, detailed=False):
    """Perform segmentation and copy number calling on normalized inputs
    """
    if not out_file:
        out_file = "%s.cns" % os.path.splitext(cnr_file)[0]
    if not utils.file_uptodate(out_file, cnr_file):
        with file_transaction(data, out_file) as tx_out_file:
            if not _cna_has_values(cnr_file):
                with open(tx_out_file, "w") as out_handle:
                    out_handle.write("chromosome\tstart\tend\tgene\tlog2\tprobes\tCN1\tCN2\tbaf\tweight\n")
            else:
                # Scale cores to avoid memory issues with segmentation
                # https://github.com/etal/cnvkit/issues/346
                if cov_interval == "genome":
                    cores = max(1, dd.get_cores(data) // 2)
                else:
                    cores = dd.get_cores(data)
                cmd = [_get_cmd(), "segment", "-p", str(cores), "-o", tx_out_file, cnr_file]
                small_vrn_files = _compatible_small_variants(data, items)
                if len(small_vrn_files) > 0 and _cna_has_values(cnr_file) and cov_interval != "genome":
                    cmd += ["--vcf", small_vrn_files[0].name, "--sample-id", small_vrn_files[0].sample]
                    if small_vrn_files[0].normal:
                        cmd += ["--normal-id", small_vrn_files[0].normal]
                resources = config_utils.get_resources("options", data["config"])
                user_options = resources.get("cnvkit_segment", [])
                cmd += [str(x) for x in user_options]
                if cov_interval == "genome" and "--threshold" not in user_options:
                    cmd += ["--threshold", "0.00001"]
                # For tumors, remove very low normalized regions, avoiding upcaptured noise
                # https://github.com/bcbio/bcbio-nextgen/issues/2171#issuecomment-348333650
                # unless we want detailed segmentation for downstream tools
                paired = vcfutils.get_paired(items)
                if paired:
                    #if detailed:
                    #    cmd += ["-m", "hmm-tumor"]
                    if "--drop-low-coverage" not in user_options:
                        cmd += ["--drop-low-coverage"]
                # preferentially use conda installed Rscript
                export_cmd = ("%s && export TMPDIR=%s && "
                              % (utils.get_R_exports(), os.path.dirname(tx_out_file)))
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
    background_file = cnvkit_background([], background_file, items, target_bed, antitarget_bed)
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
        header = next(in_handle).strip().split("\t")[1:]
        vals = map(float, next(in_handle).strip().split("\t")[1:])
    return dict(zip(header, vals))

@utils.map_wrap
@zeromq_aware_logging
def _cnvkit_fix(cnns, background_cnn, items, ckouts):
    """Normalize samples, correcting sources of bias.
    """
    return [_cnvkit_fix_base(cnns, background_cnn, items, ckouts)]

def _cnvkit_fix_base(cnns, background_cnn, items, ckouts, ext=""):
    assert len(cnns) == 2, "Expected target and antitarget CNNs: %s" % cnns
    target_cnn = [x["file"] for x in cnns if x["cnntype"] == "target"][0]
    antitarget_cnn = [x["file"] for x in cnns if x["cnntype"] == "antitarget"][0]
    dindex, data = [(i, x) for (i, x) in enumerate(items) if dd.get_sample_name(x) == cnns[0]["sample"]][0]
    out_file = ckouts[dindex]["cnr"]
    return (run_fix(target_cnn, antitarget_cnn, background_cnn, out_file, data), data)

@utils.map_wrap
@zeromq_aware_logging
def run_fix_parallel(target_cnn, antitarget_cnn, background_cnn, out_file, data):
    return [run_fix(target_cnn, antitarget_cnn, background_cnn, out_file, data)]

def run_fix(target_cnn, antitarget_cnn, background_cnn, out_file, data):
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [_get_cmd(), "fix", "-o", tx_out_file, target_cnn, antitarget_cnn, background_cnn,
                   "--sample-id", dd.get_sample_name(data)]
            do.run(_prep_cmd(cmd, tx_out_file), "CNVkit fix")
    return out_file

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

def cnvkit_background(background_cnns, out_file, items, target_bed=None, antitarget_bed=None):
    """Calculate background reference, handling flat case with no normal sample.
    """
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            cmd = [_get_cmd(), "reference", "-f", dd.get_ref_file(items[0]), "-o", tx_out_file]
            gender = _get_batch_gender(items)
            if gender:
                cmd += ["--sample-sex", gender]
            if len(background_cnns) == 0:
                assert target_bed and antitarget_bed, "Missing CNNs and target BEDs for flat background"
                cmd += ["-t", target_bed, "-a", antitarget_bed]
            else:
                cmd += background_cnns
            do.run(_prep_cmd(cmd, tx_out_file), "CNVkit background")
    return out_file


def _get_batch_gender(items):
    """Retrieve gender for a batch of items if consistent.

    Better not to specify for mixed populations, CNVkit will work
    it out
    https://github.com/bcbio/bcbio-nextgen/commit/1a0e217c8a4d3cee10fa890fb3cfd4db5034281d#r26279752
    """
    genders = set([population.get_gender(x) for x in items])
    if len(genders) == 1:
        gender = genders.pop()
        if gender != "unknown":
            return gender

def targets_w_bins(cnv_file, access_file, target_anti_fn, work_dir, data):
    """Calculate target and anti-target files with pre-determined bins.
    """
    target_file = os.path.join(work_dir, "%s-target.bed" % dd.get_sample_name(data))
    anti_file = os.path.join(work_dir, "%s-antitarget.bed" % dd.get_sample_name(data))
    if not utils.file_exists(target_file):
        target_bin, _ = target_anti_fn()
        with file_transaction(data, target_file) as tx_out_file:
            cmd = [_get_cmd(), "target", cnv_file, "--split", "-o", tx_out_file,
                   "--avg-size", str(target_bin)]
            do.run(_prep_cmd(cmd, tx_out_file), "CNVkit target")
    if not os.path.exists(anti_file):
        _, anti_bin = target_anti_fn()
        with file_transaction(data, anti_file) as tx_out_file:
            # Create access file without targets to avoid overlap
            # antitarget in cnvkit is meant to do this but appears to not always happen
            # after chromosome 1
            tx_access_file = os.path.join(os.path.dirname(tx_out_file), os.path.basename(access_file))
            pybedtools.BedTool(access_file).subtract(cnv_file).saveas(tx_access_file)
            cmd = [_get_cmd(), "antitarget", "-g", tx_access_file, cnv_file, "-o", tx_out_file,
                   "--avg-size", str(anti_bin)]
            do.run(_prep_cmd(cmd, tx_out_file), "CNVkit antitarget")
    return target_file, anti_file

def targets_from_background(back_cnn, work_dir, data):
    """Retrieve target and antitarget BEDs from background CNN file.
    """
    target_file = os.path.join(work_dir, "%s.target.bed" % dd.get_sample_name(data))
    anti_file = os.path.join(work_dir, "%s.antitarget.bed" % dd.get_sample_name(data))
    if not utils.file_exists(target_file):
        with file_transaction(data, target_file) as tx_out_file:
            out_base = tx_out_file.replace(".target.bed", "")
            cmd = [_get_cmd("reference2targets.py"), "-o", out_base, back_cnn]
            do.run(_prep_cmd(cmd, tx_out_file), "CNVkit targets from background")
            shutil.copy(out_base + ".antitarget.bed", anti_file)
    return target_file, anti_file

def _add_seg_to_output(out, data, enumerate_chroms=False):
    """Export outputs to 'seg' format compatible with IGV and GenePattern.
    """
    out_file = "%s.seg" % os.path.splitext(out["cns"])[0]
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [os.path.join(os.path.dirname(sys.executable), "cnvkit.py"), "export",
                   "seg"]
            if enumerate_chroms:
                cmd += ["--enumerate-chroms"]
            cmd += ["-o", tx_out_file, out["cns"]]
            do.run(cmd, "CNVkit export seg")
    out["seg"] = out_file
    return out

def _add_cnr_bedgraph_and_bed_to_output(out, data):
    cnr_file = out["cnr"]
    bedgraph_file = os.path.join(_sv_workdir(data), os.path.basename(cnr_file) + ".bedgraph")
    if not utils.file_exists(bedgraph_file):
        with file_transaction(data, bedgraph_file) as tx_out_file:
            cmd = "sed 1d {cnr_file} | cut -f1,2,3,5 > {tx_out_file}"
            do.run(cmd.format(**locals()), "Converting cnr to bedgraph format")
    out["cnr_bedgraph"] = bedgraph_file

    bed_file = os.path.join(_sv_workdir(data), os.path.basename(cnr_file) + ".bed")
    if not utils.file_exists(bed_file):
        with file_transaction(data, bed_file) as tx_out_file:
            cmd = "sed 1d {cnr_file} | cut -f1,2,3,4,5 > {tx_out_file}"
            do.run(cmd.format(**locals()), "Converting cnr to bed format")
    out["cnr_bed"] = bed_file
    return out

def _compatible_small_variants(data, items):
    """Retrieve small variant (SNP, indel) VCFs compatible with CNVkit.
    """
    from bcbio import heterogeneity
    VarFile = collections.namedtuple("VarFile", ["name", "sample", "normal"])
    out = []
    paired = vcfutils.get_paired(items)
    for v in heterogeneity.get_variants(data, include_germline=not paired):
        vrn_file = v["vrn_file"]
        base, ext = utils.splitext_plus(os.path.basename(vrn_file))
        if paired:
            out.append(VarFile(vrn_file, paired.tumor_name, paired.normal_name))
        else:
            out.append(VarFile(vrn_file, dd.get_sample_name(data), None))
    return out

def _add_variantcalls_to_output(out, data, items, is_somatic=False):
    """Call ploidy and convert into VCF and BED representations.
    """
    call_file = "%s-call%s" % os.path.splitext(out["cns"])
    if not utils.file_exists(call_file):
        with file_transaction(data, call_file) as tx_call_file:
            filters = ["--filter", "cn"]
            cmd = [os.path.join(os.path.dirname(sys.executable), "cnvkit.py"), "call"] + \
                  filters + \
                   ["--ploidy", str(ploidy.get_ploidy([data])),
                    "-o", tx_call_file, out["cns"]]
            small_vrn_files = _compatible_small_variants(data, items)
            if len(small_vrn_files) > 0 and _cna_has_values(out["cns"]):
                cmd += ["--vcf", small_vrn_files[0].name, "--sample-id", small_vrn_files[0].sample]
                if small_vrn_files[0].normal:
                    cmd += ["--normal-id", small_vrn_files[0].normal]
                if not is_somatic:
                    cmd += ["-m", "clonal"]
            gender = _get_batch_gender(items)
            if gender:
                cmd += ["--sample-sex", gender]
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
                do.run(cmd, "CNVkit export %s" % outformat)
    out["call_file"] = call_file
    out["vrn_bed"] = annotate.add_genes(calls["bed"], data)
    effects_vcf, _ = effects.add_to_vcf(calls["vcf"], data, "snpeff")
    out["vrn_file"] = effects_vcf or calls["vcf"]
    out["vrn_file"] = shared.annotate_with_depth(out["vrn_file"], items)
    return out

def _add_segmetrics_to_output(out, data):
    """Add metrics for measuring reliability of CNV estimates.
    """
    out_file = "%s-segmetrics.txt" % os.path.splitext(out["cns"])[0]
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [os.path.join(os.path.dirname(sys.executable), "cnvkit.py"), "segmetrics",
                   "--median", "--iqr", "--ci", "--pi",
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
            gender = _get_batch_gender([data])
            if gender:
                cmd += ["--sample-sex", gender]
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
    if len(all_sizes) > 5:
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
    else:
        thresh = 0
    larger_chroms = []
    for c in ref.file_contigs(ref_file):
        if c.size > thresh:
            larger_chroms.append(c.name)
    return larger_chroms

def _remove_haplotype_chroms(in_file, data):
    """Remove shorter haplotype chromosomes from cns/cnr files for plotting.
    """
    larger_chroms = set(_get_larger_chroms(dd.get_ref_file(data)))
    out_file = os.path.join(_sv_workdir(data), "%s-chromfilter%s" % utils.splitext_plus(os.path.basename(in_file)))
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(in_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    for line in in_handle:
                        if ((line.find("chromosome") >= 0 and line.find("start") >= 0)
                              or line.split()[0] in larger_chroms):
                            out_handle.write(line)
    return out_file

def _add_global_scatter_plot(out, data):
    out_file = os.path.join(_sv_workdir(data),
                            os.path.splitext(os.path.basename(out["cnr"]))[0] + "-scatter_global.pdf")
    if utils.file_exists(out_file):
        return out_file
    cnr = _remove_haplotype_chroms(out["cnr"], data)
    cns = _remove_haplotype_chroms(out["cns"], data)
    with file_transaction(data, out_file) as tx_out_file:
        cmd = [_get_cmd(), "scatter", "-s", cns, "-o", tx_out_file, cnr]
        do.run(_prep_cmd(cmd, tx_out_file), "CNVkit global scatter plot")
    return out_file

def _add_scatter_plot(out, data):
    out_file = os.path.join(_sv_workdir(data),
                            os.path.splitext(os.path.basename(out["cnr"]))[0] + "-scatter.pdf")
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
    out_file = os.path.join(_sv_workdir(data),
                            os.path.splitext(os.path.basename(out["cnr"]))[0] + "-diagram.pdf")
    cnr = _remove_haplotype_chroms(out["cnr"], data)
    cns = _remove_haplotype_chroms(out["cns"], data)
    if _cnx_is_empty(cnr) or _cnx_is_empty(cns):
        return None
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [_get_cmd(), "diagram", "-s", cns,
                   "-o", tx_out_file, cnr]
            gender = _get_batch_gender([data])
            if gender:
                cmd += ["--sample-sex", gender]
            do.run(_prep_cmd(cmd, tx_out_file), "CNVkit diagram plot")
    return out_file

def segment_from_cnr(cnr_file, data, out_base):
    """Provide segmentation on a cnr file, used in external PureCN integration.
    """
    cns_file = _cnvkit_segment(cnr_file, dd.get_coverage_interval(data),
                               data, [data], out_file="%s.cns" % out_base, detailed=True)
    out = _add_seg_to_output({"cns": cns_file}, data, enumerate_chroms=False)
    return out["seg"]

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
