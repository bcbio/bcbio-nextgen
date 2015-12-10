import os
import pandas as pd
import pybedtools

import numpy as np

from bcbio.utils import (file_exists, tmpfile, chdir, splitext_plus,
                         max_command_length, robust_partition_all, append_stem)
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio import broad
from bcbio.pipeline import config_utils
from bcbio.variation import vcfutils

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
    return out_file

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
                print >>out_handle, regions["header"]
                for line in in_handle:
                    idx = "".join(line.split("\t")[:2])
                    if idx not in regions:
                        print >>out_handle, "%s\t1000\t1000\t100\t100\t100\t100\t100\t100\t100\t100\t100\t100\t%s" % (line.strip(), sample)
                    else:
                        print >>out_handle, regions[idx]
    return out_file

def coverage(data):
    """
    Calculate coverage at different completeness cutoff
    for region in coverage option.
    """
    bed_file = dd.get_coverage(data)
    if not bed_file:
        return data

    work_dir = os.path.join(dd.get_work_dir(data), "report", "coverage")
    with chdir(work_dir):
        in_bam = data['work_bam']
        sample = dd.get_sample_name(data)
        logger.debug("doing coverage for %s" % sample)
        parse_file = os.path.join(sample + "_coverage.bed")
        parse_total_file = os.path.join(sample + "_cov_total.tsv")
        cores = dd.get_num_cores(data)
        if not file_exists(parse_file):
            with file_transaction(parse_file) as out_tx:
                cmd = ("sambamba depth region -F \"not unmapped\" -t {cores} -C 1000 -T 1 -T 5 -T 10 -T 20 -T 40 -T 50 -T 60 -T 70 -T 80 -T 100 -L {bed_file}  {in_bam} | sed 's/# chrom/chrom/' > {parse_file}")
                do.run(cmd.format(**locals()), "Run coverage for {}".format(sample))
        parse_file = _add_high_covered_regions(parse_file, bed_file, sample)
        _calculate_percentiles(parse_file, sample)
        data['coverage'] = os.path.abspath(parse_file)
        return data

def variants(data):
    if not "vrn_file" in  data:
        return data
    if not dd.get_coverage(data):
        return data

    in_vcf = data['vrn_file']
    work_dir = os.path.join(dd.get_work_dir(data), "report", "variants")
    with chdir(work_dir):
        in_bam = data['work_bam']
        ref_file = dd.get_ref_file(data)
        assert ref_file, "Need the reference genome fasta file."
        bed_file = dd.get_variant_regions(data)
        sample = dd.get_sample_name(data)
        in_bam = data.get("work_bam")
        cg_file = os.path.join(sample + "_with-gc.vcf.gz")
        parse_file = os.path.join(sample + "_gc-depth-parse.tsv")
        num_cores = dd.get_num_cores(data)
        broad_runner = broad.runner_from_config_safe(data["config"])
        if in_bam and broad_runner and broad_runner.has_gatk():
            if not file_exists(cg_file):
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
        return data

def priority_coverage(data):
    AVERAGE_REGION_STRING_LENGTH = 100
    bed_file = dd.get_priority_regions(data)
    if not bed_file:
        return data

    work_dir = os.path.join(dd.get_work_dir(data), "report", "coverage")
    batch_size = max_command_length() / AVERAGE_REGION_STRING_LENGTH

    sample = dd.get_sample_name(data)
    out_file = os.path.join(sample + "_priority_depth.bed")
    if file_exists(out_file):
        data['priority_coverage'] = os.path.abspath(out_file)
        return data
    with chdir(work_dir):
        in_bam = data['work_bam']
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
                    line_batch += str(line)
                if not coord_batch:
                    continue
                region_file = pybedtools.BedTool(line_batch,
                                                from_string=True).saveas().fn
                coord_string = " ".join(coord_batch)
                awk_string = r"""'BEGIN {OFS="\t"} {print $1,$2+$5,$2+$5,$4,$6"\t%s"}'""" % sample
                cmd = ("samtools view -b {in_bam} {coord_string} | "
                        "bedtools coverage -d -a {region_file} -b - | "
                        "awk {awk_string} >> {tx_out_file}")
                _silence_run(cmd.format(**locals()))
        data['priority_coverage'] = os.path.abspath(out_file)
    return data

def priority_total_coverage(data):
    """
    calculate coverage at 10 depth intervals in the priority regions
    """
    bed_file = dd.get_priority_regions(data)
    if not bed_file:
        return data
    work_dir = os.path.join(dd.get_work_dir(data), "report", "coverage")
    sample = dd.get_sample_name(data)
    out_file = os.path.join(work_dir, sample + "_priority_total_coverage.bed")
    if file_exists(out_file):
        data['priority_total_coverage'] = os.path.abspath(out_file)
        return data

    nthreads = dd.get_num_cores(data)
    in_bam = dd.get_work_bam(data)
    sambamba = config_utils.get_program("sambamba", data, default="sambamba")
    with file_transaction(out_file) as tx_out_file:
        cmd = ("{sambamba} depth region -t {nthreads} -L {bed_file} "
               "-F \"not unmapped\" "
               "-T 10 -T 20 -T 30 -T 40 -T 50 -T 60 -T 70 -T 80 -T 90 -T 100 "
               "{in_bam} -o {tx_out_file}")
        message = "Calculating coverage of {bed_file} regions in {in_bam}"
        do.run(cmd.format(**locals()), message.format(**locals()))
    data['priority_total_coverage'] = os.path.abspath(out_file)
    return data
