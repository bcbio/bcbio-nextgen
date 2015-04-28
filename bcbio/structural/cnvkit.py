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
import numpy as np
import toolz as tz

from bcbio import install, utils
from bcbio.bam import ref
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.pipeline import datadict as dd
from bcbio.variation import bedutils, vcfutils
from bcbio.provenance import do
from bcbio.structural import annotate, shared, theta, plot

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
        ckout = _add_plots_to_output(ckout, data)
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
    # Skip THetA runs until we can speed up data preparation steps
    # ckout = theta.run(ckout, paired)
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

def _get_cmd():
    return os.path.join(os.path.dirname(sys.executable), "cnvkit.py")

def _run_cnvkit_shared(data, test_bams, background_bams, access_file, work_dir,
                       background_name=None):
    """Shared functionality to run CNVkit.
    """
    ref_file = dd.get_ref_file(data)
    raw_work_dir = os.path.join(work_dir, "raw")
    out_base = os.path.splitext(os.path.basename(test_bams[0]))[0].split(".")[0]
    background_cnn = "%s_background.cnn" % (background_name if background_name else "flat")
    files = {"cnr": os.path.join(raw_work_dir, "%s.cnr" % out_base),
             "cns": os.path.join(raw_work_dir, "%s.cns" % out_base),
             "back_cnn": os.path.join(raw_work_dir, background_cnn)}
    if not utils.file_exists(files["cnr"]):
        if os.path.exists(raw_work_dir):
            shutil.rmtree(raw_work_dir)
        with tx_tmpdir(data, work_dir) as tx_work_dir:
            # pick targets, anti-targets and access files based on analysis type
            # http://cnvkit.readthedocs.org/en/latest/nonhybrid.html
            cov_interval = dd.get_coverage_interval(data)
            base_regions = dd.get_variant_regions(data)
            # For genome calls, subset to regions within 10kb of genes
            if cov_interval == "genome":
                base_regions = annotate.subset_by_genes(base_regions, data, work_dir, pad=1e4)
            raw_target_bed = bedutils.merge_overlaps(base_regions, data, out_dir=work_dir)
            target_bed = annotate.add_genes(raw_target_bed, data)

            if cov_interval == "amplicon":
                target_opts = ["--targets", target_bed, "--access", target_bed]
            elif cov_interval == "genome":
                target_opts = ["--targets", target_bed, "--access", dd.get_variant_regions(data)]
            else:
                target_opts = ["--targets", target_bed, "--access", access_file]

            cores = min(tz.get_in(["config", "algorithm", "num_cores"], data, 1),
                        len(test_bams) + len(background_bams))
            cmd = [_get_cmd(), "batch"] + \
                  test_bams + ["-n"] + background_bams + ["-f", ref_file] + \
                  target_opts + \
                  ["-d", tx_work_dir, "--split", "-p", str(cores),
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
    for ftype in ["cnr", "cns"]:
        if not os.path.exists(files[ftype]):
            raise IOError("Missing CNVkit %s file: %s" % (ftype, files[ftype]))
    return files

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
                   "freebayes", "--sample-id", dd.get_sample_name(data),
                   "--ploidy", str(dd.get_ploidy(data)),
                   "-o", tx_out_file, out["cns"]]
            gender = dd.get_gender(data)
            if gender:
                cmd += ["--gender", gender]
                if gender.lower() == "male":
                    cmd += ["--male-reference"]
            do.run(cmd, "CNVkit export FreeBayes BED cnvmap")
    out["vrn_file"] = annotate.add_genes(out_file, data)
    return out

def _add_plots_to_output(out, data):
    """Add CNVkit plots summarizing called copy number values.
    """
    out["plot"] = {"diagram": _add_diagram_plot(out, data)}
    loh_plot = _add_loh_plot(out, data)
    if loh_plot:
        out["plot"]["loh"] = loh_plot
    scatter_plot = _add_scatter_plot(out, data)
    if scatter_plot:
        out["plot"]["scatter"] = scatter_plot
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

def _add_scatter_plot(out, data):
    out_file = "%s-scatter.pdf" % os.path.splitext(out["cnr"])[0]
    priority_regions = dd.get_priority_regions(data)
    if not priority_regions:
        return None
    priority_bed = plot._prioritize_plot_regions(pybedtools.BedTool(priority_regions), data)
    if utils.file_exists(out_file):
        return out_file
    cnr = _remove_haplotype_chroms(out["cnr"], data)
    cns = _remove_haplotype_chroms(out["cns"], data)
    with file_transaction(data, out_file) as tx_out_file:
        cmd = [_get_cmd(), "scatter", "-s", cns, "-o", tx_out_file, "-l",
               priority_bed, cnr]
        do.run(cmd, "CNVkit scatter plot")
    return out_file

def _add_diagram_plot(out, data):
    out_file = "%s-diagram.pdf" % os.path.splitext(out["cnr"])[0]
    cnr = _remove_haplotype_chroms(out["cnr"], data)
    cns = _remove_haplotype_chroms(out["cns"], data)
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [_get_cmd(), "diagram", "-s", cns,
                   "-o", tx_out_file, cnr]
            gender = dd.get_gender(data)
            if gender and gender.lower() == "male":
                cmd += ["--male-reference"]
            do.run(cmd, "CNVkit diagram plot")
    return out_file

def _add_loh_plot(out, data):
    vrn_files = filter(lambda x: x is not None, [x.get("vrn_file") for x in data.get("variants", [])])
    if len(vrn_files) > 0:
        out_file = "%s-loh.pdf" % os.path.splitext(out["cnr"])[0]
        cns = _remove_haplotype_chroms(out["cns"], data)
        if not utils.file_exists(out_file):
            with file_transaction(data, out_file) as tx_out_file:
                cmd = [_get_cmd(), "loh", "-t", "-s", cns,
                       "-o", tx_out_file, vrn_files[0]]
                do.run(cmd, "CNVkit diagram plot")
        return out_file

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
    avg_size = np.median(sizes) if len(sizes) > 0 else 0
    if len(sizes) < 500 and avg_size < 10000.0:  # Default antitarget-min-size
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
