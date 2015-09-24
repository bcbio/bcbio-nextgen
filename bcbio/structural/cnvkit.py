"""Copy number detection with CNVkit with specific support for targeted sequencing.

http://cnvkit.readthedocs.org
"""
import copy
import os
import sys
import tempfile

import pybedtools
import numpy as np
import toolz as tz

from bcbio import install, utils
from bcbio.bam import ref
from bcbio.distributed.multi import run_multicore, zeromq_aware_logging
from bcbio.distributed.transaction import file_transaction
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

def _cnvkit_by_type(items, background):
    """Dispatch to specific CNVkit functionality based on input type.
    """
    if len(items + background) == 1:
        return _run_cnvkit_single(items[0])
    elif vcfutils.get_paired_phenotype(items[0]):
        return _run_cnvkit_cancer(items, background)
    else:
        return _run_cnvkit_population(items, background)

def _associate_cnvkit_out(ckouts, items):
    """Associate cnvkit output with individual items.
    """
    assert len(ckouts) == len(items)
    out = []
    for ckout, data in zip(ckouts, items):
        ckout = copy.deepcopy(ckout)
        ckout["variantcaller"] = "cnvkit"
        ckout = _add_seg_to_output(ckout, data)
        ckout = _add_gainloss_to_output(ckout, data)
        ckout = _add_segmetrics_to_output(ckout, data)
        ckout = _add_variantcalls_to_output(ckout, data)
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
    work_dir = _sv_workdir(data)
    test_bams = [data["align_bam"]]
    if background:
        background_bams = [x["align_bam"] for x in background]
        background_name = os.path.splitext(os.path.basename(background_bams[0]))[0]
    else:
        background_bams = []
        background_name = None
    ckouts = _run_cnvkit_shared([data], test_bams, background_bams, work_dir,
                               background_name=background_name)
    if not ckouts:
        return [data]
    else:
        assert len(ckouts) == 1
        return _associate_cnvkit_out(ckouts, [data])

def _run_cnvkit_cancer(items, background):
    """Run CNVkit on a tumor/normal pair.
    """
    paired = vcfutils.get_paired_bams([x["align_bam"] for x in items], items)
    work_dir = _sv_workdir(paired.tumor_data)
    ckouts = _run_cnvkit_shared([paired.tumor_data], [paired.tumor_bam], [paired.normal_bam],
                               work_dir, background_name=paired.normal_name)
    if not ckouts:
        return items
    assert len(ckouts) == 1
    tumor_data = _associate_cnvkit_out(ckouts, [paired.tumor_data])
    normal_data = [x for x in items if dd.get_sample_name(x) != paired.tumor_name]
    return tumor_data + normal_data

def _run_cnvkit_population(items, background):
    """Run CNVkit on a population of samples.

    Tries to calculate background based on case/controls, otherwise uses
    a flat background for each sample and calls independently.
    """
    assert not background
    inputs, background = shared.find_case_control(items)
    work_dir = _sv_workdir(inputs[0])
    ckouts = _run_cnvkit_shared(inputs, [x["align_bam"] for x in inputs],
                                [x["align_bam"] for x in background], work_dir,
                                background_name=dd.get_sample_name(background[0]) if len(background) > 0 else None)
    return _associate_cnvkit_out(ckouts, inputs) + background

def _get_cmd():
    return os.path.join(os.path.dirname(sys.executable), "cnvkit.py")

def _bam_to_outbase(bam_file, work_dir):
    """Convert an input BAM file into CNVkit expected output.
    """
    out_base = os.path.splitext(os.path.basename(bam_file))[0].split(".")[0]
    return os.path.join(work_dir, out_base)

def _run_cnvkit_shared(items, test_bams, background_bams, work_dir, background_name=None):
    """Shared functionality to run CNVkit, parallelizing over multiple BAM files.
    """
    raw_work_dir = utils.safe_makedir(os.path.join(work_dir, "raw"))

    background_cnn = os.path.join(raw_work_dir,
                                  "%s_background.cnn" % (background_name if background_name else "flat"))
    ckouts = []
    for test_bam in test_bams:
        out_base = _bam_to_outbase(test_bam, raw_work_dir)
        ckouts.append({"cnr": "%s.cns" % out_base,
                       "cns": "%s.cns" % out_base,
                       "back_cnn": background_cnn})
    if not utils.file_exists(ckouts[0]["cnr"]):
        data = items[0]
        cov_interval = dd.get_coverage_interval(data)
        raw_target_bed, access_bed = _get_target_access_files(cov_interval, data, work_dir)
        # bail out if we ended up with no regions
        if not utils.file_exists(raw_target_bed):
            return {}
        raw_target_bed = annotate.add_genes(raw_target_bed, data)
        parallel = {"type": "local", "cores": dd.get_cores(data), "progs": ["cnvkit"]}
        target_bed, antitarget_bed = _cnvkit_targets(raw_target_bed, access_bed, cov_interval, raw_work_dir, data)
        def _bam_to_itype(bam):
            return "background" if bam in background_bams else "evaluate"
        split_cnns = run_multicore(_cnvkit_coverage,
                                   [(bam, bed, _bam_to_itype(bam), raw_work_dir, data)
                                    for bam in test_bams + background_bams
                                    for bed in _split_bed(target_bed, data) + _split_bed(antitarget_bed, data)],
                                   data["config"], parallel)
        coverage_cnns = _merge_coverage(split_cnns, data)
        background_cnn = _cnvkit_background([x["file"] for x in coverage_cnns if x["itype"] == "background"],
                                            background_cnn, target_bed, antitarget_bed, data)
        fixed_cnrs = run_multicore(_cnvkit_fix,
                                   [(cnns, background_cnn, data) for cnns in
                                    tz.groupby("bam", [x for x in coverage_cnns
                                                       if x["itype"] == "evaluate"]).values()],
                                      data["config"], parallel)
        called_segs = run_multicore(_cnvkit_segment,
                                    [(cnr, cov_interval, data) for cnr in fixed_cnrs],
                                    data["config"], parallel)
    return ckouts

@utils.map_wrap
@zeromq_aware_logging
def _cnvkit_segment(cnr_file, cov_interval, data):
    """Perform segmentation and copy number calling on normalized inputs
    """
    out_file = "%s.cns" % os.path.splitext(cnr_file)[0]
    if not utils.file_uptodate(out_file, cnr_file):
        with file_transaction(data, out_file) as tx_out_file:
            local_sitelib = os.path.join(install.get_defaults().get("tooldir", "/usr/local"),
                                            "lib", "R", "site-library")
            cmd = [_get_cmd(), "segment", "-o", tx_out_file, "--rlibpath", local_sitelib, cnr_file]
            if cov_interval == "genome":
                cmd += ["--threshold", "0.00001"]
            # preferentially use conda installed Rscript
            export_cmd = "export PATH=%s:$PATH && " % os.path.dirname(utils.Rscript_cmd())
            do.run(export_cmd + " ".join(cmd), "CNVkit segment")
    return out_file

@utils.map_wrap
@zeromq_aware_logging
def _cnvkit_fix(cnns, background_cnn, data):
    """Normalize samples, correcting sources of bias.
    """
    assert len(cnns) == 2, "Expected target and antitarget CNNs: %s" % cnns
    target_cnn = [x["file"] for x in cnns if x["cnntype"] == "target"][0]
    antitarget_cnn = [x["file"] for x in cnns if x["cnntype"] == "antitarget"][0]
    out_file = "%scnr" % os.path.commonprefix([target_cnn, antitarget_cnn])
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [_get_cmd(), "fix", "-o", tx_out_file, target_cnn, antitarget_cnn, background_cnn]
            do.run(cmd, "CNVkit fix")
    return [out_file]

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
            do.run(cmd, "CNVkit background")
    return out_file

def _split_bed(bed_input, data):
    """Split BED file into sections for processing, allowing better multicore usage.
    """
    split_lines = 100000
    split_info = []
    base, ext = os.path.splitext(bed_input)
    base, ext2 = os.path.splitext(base)
    ext = ext2 + ext
    with open(bed_input) as in_handle:
        for cur_index, line_group in enumerate(tz.partition_all(split_lines, in_handle)):
            cur_file = "%s-%s%s" % (base, cur_index, ext)
            if not utils.file_uptodate(cur_file, bed_input):
                with file_transaction(data, cur_file) as tx_out_file:
                    with open(tx_out_file, "w") as out_handle:
                        for line in line_group:
                            out_handle.write(line)
            split_info.append({"i": cur_index, "orig": bed_input, "file": cur_file})
    if not split_info:  # empty input file
        split_info.append({"file": bed_input, "orig": bed_input})
    return split_info

def _merge_coverage(cnns, data):
    """Merge split CNN outputs into final consolidated output.
    """
    out = []
    for (out_file, _), members in tz.groupby(lambda x: (x["final_out"], x["bed_orig"]), cnns).items():
        if not utils.file_exists(out_file):
            with file_transaction(data, out_file) as tx_out_file:
                with open(tx_out_file, "w") as out_handle:
                    for i, in_file in enumerate([x["file"] for x in sorted(members, key=lambda x: x["bed_i"])]):
                        with open(in_file) as in_handle:
                            header = in_handle.readline()
                            if i == 0:
                                out_handle.write(header)
                            for line in in_handle:
                                out_handle.write(line)
        base = copy.deepcopy(members[0])
        base = tz.dissoc(base, "final_out", "bed_i", "bed_orig")
        base["file"] = out_file
        out.append(base)
    return out

@utils.map_wrap
@zeromq_aware_logging
def _cnvkit_coverage(bam_file, bed_info, input_type, work_dir, data):
    """Calculate coverage in a BED file for CNVkit.
    """
    bed_file = bed_info["file"]
    exts = {".target.bed": ("target", "targetcoverage.cnn"),
            ".antitarget.bed": ("antitarget", "antitargetcoverage.cnn")}
    assert bed_file.endswith(tuple(exts.keys())), "Unexpected BED file extension for coverage %s" % bed_file
    for orig, (cnntype, ext) in exts.items():
        if bed_file.endswith(orig):
            break
    base = _bam_to_outbase(bam_file, work_dir)
    merged_out_file = "%s.%s" % (base, ext)
    out_file = "%s-%s.%s" % (base, bed_info["i"], ext) if "i" in bed_info else merged_out_file
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [_get_cmd(), "coverage", bam_file, bed_file, "-o", tx_out_file]
            do.run(cmd, "CNVkit coverage")
    return [{"itype": input_type, "file": out_file, "bam": bam_file, "cnntype": cnntype,
             "final_out": merged_out_file, "bed_i": bed_info.get("i"), "bed_orig": bed_info["orig"]}]

def _cnvkit_targets(raw_target_bed, access_bed, cov_interval, work_dir, data):
    """Create target and antitarget regions from target and access files.
    """
    target_bed = os.path.join(work_dir, "%s.target.bed" % os.path.splitext(os.path.basename(raw_target_bed))[0])
    if not utils.file_uptodate(target_bed, raw_target_bed):
        with file_transaction(data, target_bed) as tx_out_file:
            cmd = [_get_cmd(), "target", raw_target_bed, "--split", "-o", tx_out_file]
            if cov_interval == "genome":
                cmd += ["--avg-size", "500"]
            do.run(cmd, "CNVkit target")
    antitarget_bed = os.path.join(work_dir, "%s.antitarget.bed" % os.path.splitext(os.path.basename(raw_target_bed))[0])
    if not os.path.exists(antitarget_bed):
        with file_transaction(data, antitarget_bed) as tx_out_file:
            cmd = [_get_cmd(), "antitarget", "-g", access_bed, target_bed, "-o", tx_out_file]
            do.run(cmd, "CNVkit antitarget")
    return target_bed, antitarget_bed

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
            if base_regions:
                base_regions = shared.remove_exclude_regions(base_regions, base_regions, [data])
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

def _add_variantcalls_to_output(out, data):
    """Call ploidy and convert into VCF and BED representations.
    """
    call_file = "%s-call%s" % os.path.splitext(out["cns"])
    gender = dd.get_gender(data)
    if not utils.file_exists(call_file):
        with file_transaction(data, call_file) as tx_call_file:
            cmd = [os.path.join(os.path.dirname(sys.executable), "cnvkit.py"), "call",
                   "--ploidy", str(dd.get_ploidy(data)),
                   "-o", tx_call_file, out["cns"]]
            if gender:
                cmd += ["--gender", gender]
                if gender.lower() == "male":
                    cmd += ["--male-reference"]
            do.run(cmd, "CNVkit call ploidy")
    calls = {}
    for outformat in ["bed", "vcf"]:
        out_file = "%s.%s" % (os.path.splitext(call_file)[0], outformat)
        calls[outformat] = out_file
        if not utils.file_exists(out_file):
            with file_transaction(data, out_file) as tx_out_file:
                cmd = [os.path.join(os.path.dirname(sys.executable), "cnvkit.py"), "export",
                       outformat, "--sample-id", dd.get_sample_name(data),
                       "--ploidy", str(dd.get_ploidy(data)),
                       "-o", tx_out_file, call_file]
                if gender and gender.lower() == "male":
                    cmd += ["--male-reference"]
                do.run(cmd, "CNVkit export %s" % outformat)
    out["call_file"] = call_file
    out["vrn_bed"] = annotate.add_genes(calls["bed"], data)
    out["vrn_file"] = calls["vcf"]
    return out

def _add_segmetrics_to_output(out, data):
    """Add metrics for measuring reliability of CNV estimates.
    """
    out_file = "%s-segmetrics.txt" % os.path.splitext(out["cns"])[0]
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = [os.path.join(os.path.dirname(sys.executable), "cnvkit.py"), "segmetrics",
                   "--iqr", "--ci", "--pi",
                   "-s", out["cns"], "-o", tx_out_file, out["cnr"]]
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
            do.run(cmd, "Export CNVkit calls as inputs for TheTA2")
    ckout["theta_input"] = out_file
    return ckout
