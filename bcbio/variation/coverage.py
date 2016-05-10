"""Examine and query coverage in sequencing experiments.

Provides estimates of coverage intervals based on callable regions
"""
import os
import yaml
import shutil

import pybedtools
import pandas as pd
import numpy as np

import bcbio.bed as bed
from bcbio.utils import (file_exists, chdir, max_command_length, safe_makedir,
                         robust_partition_all, append_stem, is_gzipped,
                         open_gzipsafe, symlink_plus, copy_plus)
from bcbio.bam import ref
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio import broad
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
    offtarget_thresh = 0.10  # percent of offtarget reads required to be capture (not amplification) based
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
            offtarget_stat_file = dd.get_offtarget_stats(data)
            if not offtarget_stat_file:
                offtarget_pct = 0.0
            else:
                with open(offtarget_stat_file) as in_handle:
                    stats = yaml.safe_load(in_handle)
                if float(stats["mapped"]) > 0:
                    offtarget_pct = stats["offtarget"] / float(stats["mapped"])
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

def coverage(data, out_dir):
    """
    Calculate coverage at different completeness cutoff
    for region in coverage option.
    """
    bed_file = dd.get_coverage(data)
    sambamba = config_utils.get_program("sambamba", data["config"])
    work_dir = safe_makedir(out_dir)
    if not bed_file:
        return {}
    cleaned_bed = os.path.join(work_dir, os.path.splitext(os.path.basename(bed_file))[0] + ".cleaned.bed")
    cleaned_bed = bed.decomment(bed_file, cleaned_bed)

    with chdir(work_dir):
        in_bam = dd.get_align_bam(data) or dd.get_work_bam(data)
        sample = dd.get_sample_name(data)
        logger.debug("doing coverage for %s" % sample)
        parse_file = os.path.join(sample + "_coverage.bed")
        parse_total_file = os.path.join(sample + "_cov_total.tsv")
        cores = dd.get_num_cores(data)
        if not file_exists(parse_file):
            with tx_tmpdir(data, work_dir) as tmp_dir:
                with file_transaction(parse_file) as out_tx:
                    cmd = ("{sambamba} depth region -F \"not unmapped\" -t {cores} "
                           "%s -T 1 -T 5 -T 10 -T 20 -T 40 -T 50 -T 60 -T 70 "
                           "-T 80 -T 100 -L {cleaned_bed} {in_bam} | sed 's/# "
                           "chrom/chrom/' > {out_tx}")
                    do.run(cmd.format(**locals()) % "-C 1000", "Run coverage for {}".format(sample))
        parse_file = _add_high_covered_regions(parse_file, cleaned_bed,  sample)
        parse_file = _calculate_percentiles(os.path.abspath(parse_file), sample)
    return os.path.abspath(parse_file)

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

def variants(data, out_dir):
    if "vrn_file" not in data:
        return {}
    if not dd.get_coverage(data):
        return {}

    in_vcf = data['vrn_file']
    sample = dd.get_sample_name(data)
    cg_file = os.path.join(sample + "_with-gc.vcf.gz")
    parse_file = os.path.join(sample + "_gc-depth-parse.tsv")
    qc_file = os.path.join(sample + "_bcbio_variants.txt")
    work_dir = safe_makedir(out_dir)
    with chdir(work_dir):
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
                              "-L", bed_file,
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
            if file_exists(qc_file) and file_exists(parse_file) and file_exists(cg_file):
                os.remove(cg_file)
        return qc_file

def priority_coverage(data, out_dir):
    AVERAGE_REGION_STRING_LENGTH = 100
    bed_file = dd.get_svprioritize(data)
    if not bed_file or not file_exists(bed_file):
        return data

    work_dir = safe_makedir(out_dir)
    batch_size = max_command_length() / AVERAGE_REGION_STRING_LENGTH

    sample = dd.get_sample_name(data)
    out_file = os.path.join(work_dir, sample + "_priority_depth.bed")
    if file_exists(out_file):
        # data['priority_coverage'] = os.path.abspath(out_file)
        return out_file
    with chdir(work_dir):
        in_bam = dd.get_align_bam(data) or dd.get_work_bam(data)
        logger.debug("Calculating priority coverage for %s" % sample)
        region_bed = pybedtools.BedTool(bed_file)
        with file_transaction(out_file) as tx_out_file:
            lcount = 0
            for chunk in robust_partition_all(batch_size, region_bed):
                coord_batch = []
                line_batch = ""
                for line in chunk:
                    lcount += 1
                    chrom = line.chrom
                    start = max(line.start, 0)
                    end = line.end
                    coords = "%s:%s-%s" % (chrom, start, end)
                    coord_batch.append(coords)
                    line_batch += "%s\t%s\t%s\n" % (chrom, start, end)
                if not coord_batch:
                    continue
                region_file = pybedtools.BedTool(line_batch,
                                                from_string=True).saveas().fn
                coord_string = " ".join(coord_batch)
                awk_string = r"""'BEGIN {OFS="\t"} {print $1,$2+$5,$2+$5,$4,$6"\t%s"}'""" % sample
                samtools = config_utils.get_program("samtools", data["config"])
                bedtools = config_utils.get_program("bedtools", data["config"])
                cmd = ("{samtools} view -b {in_bam} {coord_string} | "
                        "{bedtools} coverage -sorted -d -a {region_file} -b - | "
                        "awk {awk_string} >> {tx_out_file}")
                _silence_run(cmd.format(**locals()))
        # data['priority_coverage'] = os.path.abspath(out_file)
    return out_file

def priority_total_coverage(data, out_dir):
    """
    calculate coverage at 10 depth intervals in the priority regions
    """
    bed_file = dd.get_svprioritize(data)
    if not bed_file and not file_exists(bed_file):
        return {}
    work_dir = safe_makedir(out_dir)
    sample = dd.get_sample_name(data)
    out_file = os.path.join(work_dir, sample + "_priority_total_coverage.bed")
    if file_exists(out_file):
        # data['priority_total_coverage'] = os.path.abspath(out_file)
        return out_file
    nthreads = dd.get_num_cores(data)
    in_bam = dd.get_align_bam(data) or dd.get_work_bam(data)
    sambamba = config_utils.get_program("sambamba", data, default="sambamba")
    with tx_tmpdir(data, work_dir) as tmp_dir:
        cleaned_bed = os.path.join(tmp_dir, os.path.basename(bed_file))
        cleaned_bed = bed.decomment(bed_file, cleaned_bed)
        with file_transaction(out_file) as tx_out_file:
            cmd = ("{sambamba} depth region -t {nthreads} -L {cleaned_bed} "
                   "-F \"not unmapped\" "
                   "-T 10 -T 20 -T 30 -T 40 -T 50 -T 60 -T 70 -T 80 -T 90 -T 100 "
                   "{in_bam} -o {tx_out_file}")
            message = "Calculating coverage of {bed_file} regions in {in_bam}"
            do.run(cmd.format(**locals()), message.format(**locals()))
    # data['priority_total_coverage'] = os.path.abspath(out_file)
    return out_file
