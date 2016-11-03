"""Examine and query coverage in sequencing experiments.

Provides estimates of coverage intervals based on callable regions
"""
import itertools
import os
import shutil
import yaml

import pybedtools
import pandas as pd
import numpy as np
import pysam

from bcbio.variation.bedutils import clean_file
from bcbio.utils import (file_exists, chdir, safe_makedir,
                         append_stem, is_gzipped, remove_plus,
                         open_gzipsafe, copy_plus, splitext_plus)
from bcbio.bam import ref
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio import bam, broad, utils
from bcbio.bam import sambamba
from bcbio.pipeline import config_utils
from bcbio.variation import vcfutils

def assign_interval(data):
    """Identify coverage based on percent of genome covered and relation to targets.

    Classifies coverage into 3 categories:
      - genome: Full genome coverage
      - regional: Regional coverage, like exome capture, with off-target reads
      - amplicon: Amplication based regional coverage without off-target reads
    """
    genome_cov_thresh = 0.40  # percent of genome covered for whole genome analysis
    offtarget_thresh = 0.05  # percent of offtarget reads required to be capture (not amplification) based
    if not dd.get_coverage_interval(data):
        vrs = dd.get_variant_regions(data)
        callable_file = dd.get_sample_callable(data)
        if vrs:
            seq_size = pybedtools.BedTool(vrs).total_coverage()
        else:
            seq_size = pybedtools.BedTool(callable_file).total_coverage()
        total_size = sum([c.size for c in ref.file_contigs(dd.get_ref_file(data), data["config"])])
        genome_cov_pct = seq_size / float(total_size)
        if genome_cov_pct > genome_cov_thresh:
            cov_interval = "genome"
            offtarget_pct = 0.0
        else:
            offtarget_stats_file = dd.get_offtarget_stats(data)
            if not offtarget_stats_file:
                offtarget_pct = 0.0
            else:
                with open(offtarget_stats_file) as in_handle:
                    stats = yaml.safe_load(in_handle)
                if stats.get("offtarget") and stats["mapped_unique"]:
                    offtarget_pct = float(stats["offtarget"]) / stats["mapped_unique"]
                else:
                    offtarget_pct = 0.0
            if offtarget_pct > offtarget_thresh:
                cov_interval = "regional"
            else:
                cov_interval = "amplicon"
        logger.info("%s: Assigned coverage as '%s' with %.1f%% genome coverage and %.1f%% offtarget coverage"
                    % (dd.get_sample_name(data), cov_interval, genome_cov_pct * 100.0, offtarget_pct * 100.0))
        data["config"]["algorithm"]["coverage_interval"] = cov_interval
    return data

def calculate(bam_file, data):
    """Calculate coverage in parallel using samtools depth through goleft.

    samtools depth removes duplicates and secondary reads from the counts:
    if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
    """
    params = {"window_size": 5000, "parallel_window_size": 1e5, "min": dd.get_coverage_depth_min(data),
              "high_multiplier": 20}
    prefix = os.path.join(
        utils.safe_makedir(os.path.join(dd.get_work_dir(data), "align", dd.get_sample_name(data))),
        "%s-coverage" % (dd.get_sample_name(data)))
    out_file = prefix + ".depth.bed"
    callable_file = prefix + ".callable.bed"
    variant_regions = dd.get_variant_regions_merged(data)
    variant_regions_avg_cov = get_average_coverage(data, bam_file, variant_regions,
                                                   "variant_regions", file_prefix=prefix)
    if not utils.file_uptodate(out_file, bam_file):
        ref_file = dd.get_ref_file(data)
        cmd = ["goleft", "depth", "--windowsize", str(params["window_size"]), "--q", "1",
               "--mincov", str(params["min"]), "--reference", ref_file,
               "--processes", str(dd.get_num_cores(data)), "--stats", "--ordered"]
        if variant_regions:
            window_file = "%s-tocalculate-windows.bed" % utils.splitext_plus(out_file)[0]
            if not utils.file_uptodate(window_file, bam_file):
                with file_transaction(data, window_file) as tx_out_file:
                    pybedtools.BedTool().window_maker(w=params["parallel_window_size"],
                                                      b=pybedtools.BedTool(variant_regions)).saveas(tx_out_file)
            cmd += ["--bed", window_file]
        max_depth = _get_max_depth(variant_regions_avg_cov, params, data)
        if max_depth:
            cmd += ["--maxmeandepth", str(int(max_depth))]
        with file_transaction(data, out_file) as tx_out_file:
            with utils.chdir(os.path.dirname(tx_out_file)):
                tx_callable_file = tx_out_file.replace(".depth.bed", ".callable.bed")
                prefix = tx_out_file.replace(".depth.bed", "")
                cmd += ["--prefix", prefix, bam_file]
                do.run(cmd, "Calculate coverage: %s" % dd.get_sample_name(data))
                shutil.move(tx_callable_file, callable_file)
    return out_file, callable_file, _extract_highdepth(callable_file, data), variant_regions_avg_cov

def _extract_highdepth(callable_file, data):
    out_file = callable_file.replace(".callable.bed", ".highdepth.bed")
    if not utils.file_uptodate(out_file, callable_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(callable_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    for line in in_handle:
                        parts = line.strip().split("\t")
                        if "EXCESSIVE_COVERAGE" in parts:
                            out_handle.write("\t".join(parts[:3] + ["highdepth"]) + "\n")
    return out_file

def _get_max_depth(average_coverage, params, data):
    """Calculate maximum depth based on a rough multiplier of average coverage.
    """
    if dd.get_coverage_interval(data) == "genome":
        avg_cov = max(30.0, average_coverage)
        return avg_cov * params["high_multiplier"]

def get_average_coverage(data, bam_file, bed_file=None, target_name="genome", file_prefix=None):
    logger.debug("Calculation average coverage of " + bam_file +
                 " on " + target_name + ((" " + bed_file) if bed_file else ""))
    file_prefix = file_prefix or os.path.join(
        utils.safe_makedir(os.path.join(dd.get_work_dir(data), "align", dd.get_sample_name(data))),
        "%s-coverage" % (dd.get_sample_name(data)))
    cache_file = file_prefix + "-" + target_name + "-stats.yaml"
    if utils.file_uptodate(cache_file, bam_file):
        with open(cache_file) as in_handle:
            stats = yaml.safe_load(in_handle)
        return stats["avg_coverage"]
    if bed_file:
        avg_cov = _average_target_coverage(data, bed_file, bam_file, target_name=target_name)
    else:
        avg_cov = _average_genome_coverage(data, bam_file)
    stats = {"avg_coverage": avg_cov}
    with open(cache_file, "w") as out_handle:
        yaml.safe_dump(stats, out_handle, default_flow_style=False, allow_unicode=False)
    return avg_cov

def _average_genome_coverage(data, bam_file):
    total = sum([c.size for c in ref.file_contigs(dd.get_ref_file(data), data["config"])])
    read_counts = sambamba.number_of_mapped_reads(data, bam_file, keep_dups=False)
    with pysam.Samfile(bam_file, "rb") as pysam_bam:
        read_size = np.median(list(itertools.islice((a.query_length for a in pysam_bam.fetch()), 1e5)))
    avg_cov = float(read_counts * read_size) / total
    return avg_cov

def _average_target_coverage(data, bed_file, bam_file, target_name):
    sambamba_depth_file = regions_coverage(data, bed_file, bam_file, target_name)
    avg_covs = []
    mean_cov_col = None
    total_len = 0
    with open(sambamba_depth_file) as fh:
        for line in fh:
            if line.startswith('#'):
                mean_cov_col = line.split('\t').index('meanCoverage')
                continue
            line_tokens = line.replace('\n', '').split()
            start, end = map(int, line_tokens[1:3])
            size = end - start
            avg_covs.append(float(line_tokens[mean_cov_col]) * size)
            total_len += size
    avg_cov = sum(avg_covs) / total_len if total_len > 0 else 0
    return avg_cov

def decorate_problem_regions(query_bed, problem_bed_dir):
    """
    decorate query_bed with percentage covered by BED files of regions specified
    in the problem_bed_dir
    """
    if is_gzipped(query_bed):
        stem, _ = os.path.splitext(query_bed)
        stem, ext = os.path.splitext(stem)
    else:
        stem, ext = os.path.splitext(query_bed)
    out_file = stem + ".problem_annotated" + ext + ".gz"
    if file_exists(out_file):
        return out_file
    bed_files = _find_bed_files(problem_bed_dir)
    bed_file_string = " ".join(bed_files)
    names = [os.path.splitext(os.path.basename(x))[0] for x in bed_files]
    names_string = " ".join(names)
    with open_gzipsafe(query_bed) as in_handle:
        header = map(str, in_handle.next().strip().split())
    header = "\t".join(header + names)
    cmd = ("bedtools annotate -i {query_bed} -files {bed_file_string} "
           "-names {names_string} | sed -s 's/^#.*$/{header}/' | bgzip -c > {tx_out_file}")
    with file_transaction(out_file) as tx_out_file:
        message = "Annotate %s with problem regions." % query_bed
        do.run(cmd.format(**locals()), message)
    return out_file

def _find_bed_files(path):
    """
    recursively walk directories to find all of the BED files in the
    problem regions directory
    """
    bed_files = []
    for dirpath, subdirs, files in os.walk(path):
        for x in files:
            if x.endswith(".bed") or x.endswith(".bed.gz"):
                bed_files.append(os.path.join(dirpath, x))
    return bed_files

def _silence_run(cmd):
    do._do_run(cmd, False)

def checkpoint(stem):
    def check_file(f):
        def wrapper(*args, **kwargs):
            out_file = append_stem(args[0], stem)
            if file_exists(out_file):
                logger.debug("Skipping %s" % out_file)
                return out_file
            return f(*args, **kwargs)
        return wrapper
    return check_file

@checkpoint("_summary")
def _calculate_percentiles(in_file, sample):
    """
    Parse pct bases per region to summarize it in
    7 different pct of regions points with pct bases covered
    higher than a completeness cutoff (5, 10, 20, 50 ...)
    """
    has_data = False
    with open(in_file) as in_handle:
        for i, line in enumerate(in_handle):
            if i > 0:
                has_data = True
                break
    if not has_data:
        return in_file
    out_file = append_stem(in_file, "_summary")
    out_total_file = append_stem(in_file, "_total_summary")
    dt = pd.read_csv(in_file, sep="\t", index_col=False)
    pct = dict()
    pct_bases = dict()
    size = np.array(dt["chromEnd"]) - np.array(dt["chromStart"])
    for cutoff in [h for h in list(dt) if h.startswith("percentage")]:
        a = np.array(dt[cutoff])
        for p_point in [0.01, 10, 25, 50, 75, 90, 99.9]:
            q = np.percentile(a, p_point)
            pct[(cutoff, p_point)] = q
        pct_bases[cutoff] = sum(size * a)/float(sum(size))

    with file_transaction(out_total_file) as tx_file:
        with open(tx_file, 'w') as out_handle:
            print >>out_handle, "cutoff_reads\tbases_pct\tsample"
            for k in pct_bases:
                print >>out_handle, "\t".join(map(str, [k, pct_bases[k], sample]))
    with file_transaction(out_file) as tx_file:
        with open(tx_file, 'w') as out_handle:
            print >>out_handle, "cutoff_reads\tregion_pct\tbases_pct\tsample"
            for k in pct:
                print >>out_handle, "\t".join(map(str, [k[0], k[1], pct[k], sample]))
    # To move metrics to multiqc, will remove older files
    # when bcbreport accepts these one, to avoid errors
    # while porting everything to multiqc
    # These files will be copied to final
    out_file_fixed = os.path.join(os.path.dirname(out_file), "%s_bcbio_coverage.txt" % sample)
    out_total_fixed = os.path.join(os.path.dirname(out_file), "%s_bcbio_coverage_avg.txt" % sample)
    copy_plus(out_file, out_file_fixed)
    copy_plus(out_total_file, out_total_fixed)
    return out_file_fixed

def _read_regions(fn):
    """
    Save in a dict the position of regions with
    the information of the coverage stats.
    """
    regions = {}
    with open(fn) as in_handle:
        for line in in_handle:
            if line.startswith("chrom"):
                regions["header"] = line.strip()
                continue
            idx = "".join(line.split("\t")[:2])
            regions[idx] = line.strip()
    return regions

@checkpoint("_fixed")
def _add_high_covered_regions(in_file, bed_file, sample):
    """
    Add regions with higher coverage than the limit
    as fully covered.
    """
    out_file = append_stem(in_file, "_fixed")
    regions = _read_regions(in_file)
    with file_transaction(out_file) as out_tx:
        with open(bed_file) as in_handle:
            with open(out_tx, 'w') as out_handle:
                if "header" in regions:
                    print >>out_handle, regions["header"]
                for line in in_handle:
                    idx = "".join(line.split("\t")[:2])
                    if idx not in regions:
                        print >>out_handle, "%s\t1000\t1000\t100\t100\t100\t100\t100\t100\t100\t100\t100\t100\t%s" % (line.strip(), sample)
                    else:
                        print >>out_handle, regions[idx]
    return out_file

def _summary_variants(in_file, out_file):
    """Parse GC and depth variant file
       to be ready for multiqc.
    """
    dt = pd.read_csv(in_file, sep="\t", index_col=False,
                     dtype={"CG": np.float64, "depth": np.float64}, na_values=["."]).dropna()
    row = list()
    with file_transaction(out_file) as out_tx:
        cg = dt["CG"]
        d = dt["depth"]
        for p_point in [0.01, 10, 25, 50, 75, 90, 99.9, 100]:
            if len(cg) > 0:
                q_cg = np.percentile(cg, p_point)
            else:
                q_cg = 0
            if len(d) > 0:
                q_d = np.percentile(d, p_point)
            else:
                q_d = 0
            row.append([p_point, q_d, q_cg])
        pd.DataFrame(row).to_csv(out_tx, header=["pct_variants", "depth", "cg"], index=False, sep="\t")

def _read_bcffile(out_file):
    out = {}
    with open(out_file) as in_handle:
        for line in in_handle:
            if line.startswith("SN") and line.find("records") > -1:
                out["Variations (total)"] = line.split()[-1]
            elif line.startswith("SN") and line.find("SNPs") > -1:
                out["Variations (SNPs)"] = line.split()[-1]
            elif line.startswith("SN") and line.find("indels") > -1:
                out["Variations (indels)"] = line.split()[-1]
            elif line.startswith("TSTV"):
                out["Variations (ts/tv)"] = line.split()[4]
            elif line.startswith("PSC"):
                # out["Variations (homozygous)"] = line.split()[3]
                out["Variations (alt homozygous)"] = line.split()[4]
                out["Variations (heterozygous)"] = line.split()[5]
    return out

def _get_variant_callers(data):
    """Use first caller if ensemble is not active"""
    callers = dd.get_variantcaller(data)
    if not callers:
        return None
    if isinstance(callers, basestring):
        callers = [callers]
    active_callers = [c.get("variantcaller") for c in data.get("variants", [{}])]
    active_vcf = [c.get("vrn_file") for c in data.get("variants", [{}])]
    active_germline = [c.get("germline") for c in data.get("variants", [{}])]
    vcf = dict(zip(active_callers, active_vcf))
    germline = dict(zip(active_callers, active_germline))
    if "ensemble" in active_callers:
        vcf_fn = vcf["ensemble"]
    else:
        vcf_fn = vcf[callers[0]]
    if not vcf_fn:
        vcf_fn = germline[callers[0]]
    return vcf_fn

def _run_bcftools(data, out_dir):
    """Get variants stats"""
    vcf_file = _get_variant_callers(data)
    opts = "-f PASS"
    out = {}
    if vcf_file:
        name = dd.get_sample_name(data)
        stem = os.path.join(out_dir, os.path.basename(splitext_plus(vcf_file)[0]))
        out_file = "%s-%s-bcfstats.tsv" % (stem, name)
        bcftools = config_utils.get_program("bcftools", data["config"])
        if not file_exists(out_file):
            cmd = ("{bcftools} stats -s {name} {opts} {vcf_file} > {out_file}")
            do.run(cmd.format(**locals()), "basic vcf stats %s" % dd.get_sample_name(data))
        out[name] = _read_bcffile(out_file)
    return out

def variants(data, out_dir):
    """Variants QC metrics"""
    if not "variants" in data:
        return None
    work_dir = safe_makedir(out_dir)
    sample = dd.get_sample_name(data)
    bcfstats = _run_bcftools(data, work_dir)
    bed_file = dd.get_coverage(data)
    bcf_out = os.path.join(sample + "_bcbio_variants_stats.txt")
    cg_file = os.path.join(sample + "_with-gc.vcf.gz")
    parse_file = os.path.join(sample + "_gc-depth-parse.tsv")
    qc_file = os.path.join(sample + "_bcbio_variants.txt")
    with chdir(work_dir):
        if not file_exists(bcf_out):
            with open(bcf_out, "w") as out_handle:
                yaml.safe_dump(bcfstats, out_handle, default_flow_style=False, allow_unicode=False)
        if "vrn_file" not in data or not bed_file:
            return None

        in_vcf = data['vrn_file']
        cleaned_bed = clean_file(bed_file, data)
        if file_exists(qc_file):
            return qc_file
        in_bam = dd.get_align_bam(data) or dd.get_work_bam(data)
        ref_file = dd.get_ref_file(data)
        assert ref_file, "Need the reference genome fasta file."
        bed_file = dd.get_variant_regions(data)
        in_bam = dd.get_align_bam(data) or dd.get_work_bam(data)
        num_cores = dd.get_num_cores(data)
        broad_runner = broad.runner_from_config_safe(data["config"])
        if in_bam and broad_runner and broad_runner.has_gatk():
            if not file_exists(parse_file):
                with file_transaction(cg_file) as tx_out:
                    params = ["-T", "VariantAnnotator",
                              "-R", ref_file,
                              "-L", cleaned_bed,
                              "-I", in_bam,
                              "-A", "GCContent",
                              "-A", "Coverage",
                              "--variant", in_vcf,
                              "--out", tx_out]
                    broad_runner.run_gatk(params)
                cg_file = vcfutils.bgzip_and_index(cg_file, data["config"])

            if not file_exists(parse_file):
                with file_transaction(parse_file) as out_tx:
                    with open(out_tx, 'w') as out_handle:
                        print >>out_handle, "CG\tdepth\tsample"
                    cmd = ("bcftools query -s {sample} -f '[%GC][\\t%DP][\\t%SAMPLE]\\n' -R "
                            "{bed_file} {cg_file} >> {out_tx}")
                    do.run(cmd.format(**locals()),
                            "Calculating GC content and depth for %s" % in_vcf)
                    logger.debug('parsing coverage: %s' % sample)
            if not file_exists(qc_file):
                # This files will be copied to final
                _summary_variants(parse_file, qc_file)
            if file_exists(qc_file) and file_exists(parse_file):
                remove_plus(cg_file)

def regions_coverage(data, bed_file, bam_file, target_name):
    work_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "coverage", dd.get_sample_name(data)))
    out_file = os.path.join(work_dir, target_name + "_regions_depth.bed")
    if utils.file_uptodate(out_file, bam_file) and utils.file_uptodate(out_file, bed_file):
        return out_file
    with file_transaction(out_file) as tx_out_file:
        cmdl = sambamba.make_command(data, "depth region", bam_file, bed_file) + " -o " + tx_out_file
        message = "Calculating regions coverage of {target_name} in {bam_file}"
        do.run(cmdl, message.format(**locals()))
    return out_file

def priority_coverage(data, out_dir):
    from bcbio.structural import prioritize
    bed_file = dd.get_svprioritize(data)
    if not bed_file or not file_exists(bed_file) or prioritize.is_gene_list(bed_file):
        return data

    work_dir = safe_makedir(out_dir)
    sample = dd.get_sample_name(data)
    cleaned_bed = clean_file(bed_file, data, prefix="cov-", simple=True)
    out_file = os.path.join(work_dir, sample + "_priority_depth.bed")
    in_bam = dd.get_align_bam(data) or dd.get_work_bam(data)
    if utils.file_uptodate(out_file, cleaned_bed) and utils.file_uptodate(out_file, in_bam):
        return out_file
    with file_transaction(out_file) as tx_out_file:
        cmdl = sambamba.make_command(data, "depth base", in_bam, cleaned_bed)
        parse_cmd = "awk '{print $1\"\t\"$2\"\t\"$2\"\t\"$3\"\t\"$10}' | sed '1d'"
        cmdl += " | {parse_cmd} > {tx_out_file}"
        message = "Calculating base coverage of {bed_file} in {in_bam}"
        do.run(cmdl.format(**locals()), message.format(**locals()))
    return out_file

def priority_total_coverage(data, out_dir):
    """
    calculate coverage at 10 depth intervals in the priority regions
    """
    from bcbio.structural import prioritize
    bed_file = dd.get_svprioritize(data)
    if not bed_file and not file_exists(bed_file) or prioritize.is_gene_list(bed_file):
        return {}
    in_bam = dd.get_align_bam(data) or dd.get_work_bam(data)
    cleaned_bed = clean_file(bed_file, data)
    work_dir = safe_makedir(out_dir)
    sample = dd.get_sample_name(data)
    out_file = os.path.join(work_dir, sample + "_priority_total_coverage.bed")
    if utils.file_uptodate(out_file, cleaned_bed) and utils.file_uptodate(out_file, in_bam):
        return out_file
    cmdl = sambamba.make_command(data, "depth region", in_bam, cleaned_bed,
                                 depth_thresholds=[10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    with file_transaction(out_file) as tx_out_file:
        message = "Calculating region coverage of {bed_file} in {in_bam}"
        do.run(cmdl + " -o " + tx_out_file, message.format(**locals()))
    logger.debug("Saved svprioritize coverage into " + out_file)
    return out_file

def coverage_region_detailed_stats(data, out_dir):
    """
    Calculate coverage at different completeness cutoff
    for region in coverage option.
    """
    bed_file = dd.get_coverage(data)
    if not bed_file:
        return None
    work_dir = safe_makedir(out_dir)
    cleaned_bed = clean_file(bed_file, data, prefix="cov-", simple=True)

    with chdir(work_dir):
        in_bam = dd.get_align_bam(data) or dd.get_work_bam(data)
        sample = dd.get_sample_name(data)
        logger.debug("doing coverage for %s" % sample)
        parse_total_file = os.path.join(sample + "_cov_total.tsv")
        parse_file = os.path.join(sample + "_coverage.bed")
        if utils.file_uptodate(parse_file, cleaned_bed) and utils.file_uptodate(parse_file, in_bam):
            pass
        else:
            with file_transaction(parse_file) as out_tx:
                cmdl = sambamba.make_command(data, "depth region", in_bam, cleaned_bed,
                                             depth_thresholds=[1, 5, 10, 20, 40, 50, 60, 70, 80, 100],
                                             max_cov=1000)
                cmdl += " | sed 's/# chrom/chrom/' > " + out_tx
                do.run(cmdl, "Run coverage regional analysis for {}".format(sample))
        parse_file = _add_high_covered_regions(parse_file, cleaned_bed, sample)
        parse_file = _calculate_percentiles(os.path.abspath(parse_file), sample)
    return os.path.abspath(parse_file)
