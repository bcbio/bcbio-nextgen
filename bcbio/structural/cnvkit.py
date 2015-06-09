"""Copy number detection with CNVkit with specific support for targeted sequencing.

http://cnvkit.readthedocs.org
"""
import os
import shutil
import sys
import tempfile

import pybedtools
import numpy as np
import toolz as tz

from bcbio import install, utils
from bcbio.bam import ref
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.heterogeneity import chromhacks
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.variation import bedutils, vcfutils
from bcbio.provenance import do
from bcbio.structural import annotate, shared, regions, plot

def run(items, background=None):
    """Detect copy number variations from batched set of samples using CNVkit.
    """
    if not background: background = []
    return _cnvkit_by_type(items, background)

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "cnvkit"))

def export_theta(ckout, data):
    """Provide updated set of data with export information for TheTA2 input.
    """
    cns_file = chromhacks.bed_to_standardonly(ckout["cns"], data, headers="chromosome")
    cnr_file = chromhacks.bed_to_standardonly(ckout["cnr"], data, headers="chromosome")
    out_file = "%s-theta.input" % utils.splitext_plus(cns_file)[0]
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [_get_cmd(), "export", "theta", cns_file, cnr_file, "-o", tx_out_file]
            do.run(cmd, "Export CNVkit calls as inputs for TheTA2")
    # ckout["theta_input"] = _subset_theta_to_calls(out_file, ckout, data)
    ckout["theta_input"] = out_file
    return ckout

def _subset_theta_to_calls(in_file, ckout, data):
    """Subset CNVkit regions to provide additional signal for THetA.

    THetA has default assumptions about lengths of calls and finding
    useful signal in longer regions. We adjust for this by subsetting
    calls to a range around the most useful signal.
    """
    tn_ratio = 0.9
    keep_background = False
    out_file = "%s-cnvsize%s" % utils.splitext_plus(in_file)
    if not utils.file_uptodate(out_file, in_file):
        call_sizes = []
        calls = set([])
        with open(ckout["vrn_file"]) as in_handle:
            for line in in_handle:
                chrom, start, end, _, count = line.split()[:5]
                if max([int(x) for x in count.split(",")]) < 6:
                    call_sizes.append((int(end) - int(start)))
                    calls.add((chrom, start, end))
        keep_min = np.percentile(call_sizes, 10)
        keep_max = np.percentile(call_sizes, 90)
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                # Pull out calls that have tumor/normal differences
                tn_count = 0
                with open(in_file) as in_handle:
                    for line in in_handle:
                        if line.startswith("#"):
                            out_handle.write(line)
                        else:
                            key = tuple(line.split()[1:4])
                            sizes = [float(x) for x in line.split()[4:6]]
                            size = int(key[2]) - int(key[1])
                            if size >= keep_min and size <= keep_max:
                                if (min(sizes) / max(sizes)) < tn_ratio:
                                    tn_count += 1
                                    out_handle.write(line)
                if keep_background:
                    # Pull out equal number of background calls
                    no_tn_count = 0
                    with open(in_file) as in_handle:
                        for line in in_handle:
                            if not line.startswith("#"):
                                key = tuple(line.split()[1:4])
                                sizes = [float(x) for x in line.split()[4:6]]
                                size = int(key[2]) - int(key[1])
                                if size >= keep_min and size <= keep_max:
                                    if no_tn_count < tn_count and (min(sizes) / max(sizes)) > tn_ratio:
                                        no_tn_count += 1
                                        out_handle.write(line)
    return out_file

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
        #  ckout = _add_coverage_bedgraph_to_output(ckout, data)
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
    work_dir = _sv_workdir(data)
    test_bams = [data["align_bam"]]
    if background:
        background_bams = [x["align_bam"] for x in background]
        background_name = os.path.splitext(os.path.basename(background_bams[0]))[0]
    else:
        background_bams = []
        background_name = None
    ckout = _run_cnvkit_shared(data, test_bams, background_bams, work_dir,
                               background_name=background_name)
    if not ckout:
        return [data]
    else:
        return _associate_cnvkit_out(ckout, [data])

def _run_cnvkit_cancer(items, background):
    """Run CNVkit on a tumor/normal pair.
    """
    paired = vcfutils.get_paired_bams([x["align_bam"] for x in items], items)
    work_dir = _sv_workdir(paired.tumor_data)
    ckout = _run_cnvkit_shared(paired.tumor_data, [paired.tumor_bam], [paired.normal_bam],
                               work_dir, background_name=paired.normal_name)
    if not ckout:
        return items

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
    return [_run_cnvkit_single(data, background)[0] for data in inputs] + \
           [_run_cnvkit_single(data, inputs)[0] for data in background]

def _get_cmd():
    return os.path.join(os.path.dirname(sys.executable), "cnvkit.py")

def _run_cnvkit_shared(data, test_bams, background_bams, work_dir, background_name=None):
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
            cov_interval = dd.get_coverage_interval(data)
            raw_target_bed, access_bed = _get_target_access_files(cov_interval, data, work_dir)
            # bail out if we ended up with no regions
            if not utils.file_exists(raw_target_bed):
                return {}
            target_bed = annotate.add_genes(raw_target_bed, data)

            cores = min(tz.get_in(["config", "algorithm", "num_cores"], data, 1),
                        len(test_bams) + len(background_bams))
            cmd = [_get_cmd(), "batch"] + \
                  test_bams + ["-n"] + background_bams + ["-f", ref_file] + \
                  ["--targets", target_bed, "--access", access_bed] + \
                  ["-d", tx_work_dir, "--split", "-p", str(cores),
                   "--output-reference", os.path.join(tx_work_dir, background_cnn)]
            if cov_interval not in ["amplicon", "genome"]:
                at_avg, at_min, t_avg = _get_antitarget_size(access_bed, target_bed)
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

def _get_target_access_files(cov_interval, data, work_dir):
    """Retrieve target and access files based on the type of data to process.

    pick targets, anti-targets and access files based on analysis type
    http://cnvkit.readthedocs.org/en/latest/nonhybrid.html
    """
    base_regions = regions.get_sv_bed(data)
    # if we don't have a configured BED or regions to use for SV caling
    if not base_regions:
        # For genome calls, subset to regions within 10kb of genes
        if cov_interval == "genome":
            base_regions = regions.get_sv_bed(data, "transcripts1e4", work_dir)
        # Finally, default to the defined variant regions
        if not base_regions:
            base_regions = dd.get_variant_regions(data)

    target_bed = bedutils.merge_overlaps(base_regions, data, out_dir=work_dir)
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
