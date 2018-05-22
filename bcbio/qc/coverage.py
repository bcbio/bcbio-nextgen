"""Coverage based QC calculations.
"""
import glob
import os
import subprocess

from bcbio.bam import ref, readstats, utils
from bcbio.distributed import transaction
from bcbio.heterogeneity import chromhacks
import bcbio.pipeline.datadict as dd
from bcbio.provenance import do
from bcbio.variation import coverage as cov
from bcbio.variation import bedutils

def run(bam_file, data, out_dir):
    """Run coverage QC analysis
    """
    out = dict()

    out_dir = utils.safe_makedir(out_dir)
    if dd.get_coverage(data) and dd.get_coverage(data) not in ["None"]:
        merged_bed_file = bedutils.clean_file(dd.get_coverage_merged(data), data, prefix="cov-", simple=True)
        target_name = "coverage"
    elif dd.get_coverage_interval(data) != "genome":
        merged_bed_file = dd.get_variant_regions_merged(data) or dd.get_sample_callable(data)
        target_name = "variant_regions"
    else:
        merged_bed_file = None
        target_name = "genome"

    avg_depth = cov.get_average_coverage(target_name, merged_bed_file, data)
    if target_name == "coverage":
        out_files = cov.coverage_region_detailed_stats(target_name, merged_bed_file, data, out_dir)
    else:
        out_files = []

    out['Avg_coverage'] = avg_depth

    samtools_stats_dir = os.path.join(out_dir, os.path.pardir, 'samtools')
    from bcbio.qc import samtools
    samtools_stats = samtools.run(bam_file, data, samtools_stats_dir)["metrics"]

    out["Total_reads"] = total_reads = int(samtools_stats["Total_reads"])
    out["Mapped_reads"] = mapped = int(samtools_stats["Mapped_reads"])
    out["Mapped_paired_reads"] = int(samtools_stats["Mapped_paired_reads"])
    out['Duplicates'] = dups = int(samtools_stats["Duplicates"])

    if total_reads:
        out["Mapped_reads_pct"] = 100.0 * mapped / total_reads
    if mapped:
        out['Duplicates_pct'] = 100.0 * dups / mapped

    if dd.get_coverage_interval(data) == "genome":
        mapped_unique = mapped - dups
    else:
        mapped_unique = readstats.number_of_mapped_reads(data, bam_file, keep_dups=False)
    out['Mapped_unique_reads'] = mapped_unique

    if merged_bed_file:
        ontarget = readstats.number_of_mapped_reads(
            data, bam_file, keep_dups=False, bed_file=merged_bed_file, target_name=target_name)
        out["Ontarget_unique_reads"] = ontarget
        if mapped_unique:
            out["Ontarget_pct"] = 100.0 * ontarget / mapped_unique
            out['Offtarget_pct'] = 100.0 * (mapped_unique - ontarget) / mapped_unique
            if dd.get_coverage_interval(data) != "genome":
                # Skip padded calculation for WGS even if the "coverage" file is specified
                # the padded statistic makes only sense for exomes and panels
                padded_bed_file = bedutils.get_padded_bed_file(out_dir, merged_bed_file, 200, data)
                ontarget_padded = readstats.number_of_mapped_reads(
                    data, bam_file, keep_dups=False, bed_file=padded_bed_file, target_name=target_name + "_padded")
                out["Ontarget_padded_pct"] = 100.0 * ontarget_padded / mapped_unique
        if total_reads:
            out['Usable_pct'] = 100.0 * ontarget / total_reads

    indexcov_files = _goleft_indexcov(bam_file, data, out_dir)
    out_files += [x for x in indexcov_files if x and utils.file_exists(x)]
    out = {"metrics": out}
    if len(out_files) > 0:
        out["base"] = out_files[0]
        out["secondary"] = out_files[1:]
    return out

def _goleft_indexcov(bam_file, data, out_dir):
    """Use goleft indexcov to estimate coverage distributions using BAM index.

    Only used for whole genome runs as captures typically don't have enough data
    to be useful for index-only summaries.
    """
    if not dd.get_coverage_interval(data) == "genome":
        return []
    out_dir = utils.safe_makedir(os.path.join(out_dir, "indexcov"))
    out_files = [os.path.join(out_dir, "%s-indexcov.%s" % (dd.get_sample_name(data), ext))
                 for ext in ["roc", "ped", "bed.gz"]]
    if not utils.file_uptodate(out_files[-1], bam_file):
        with transaction.tx_tmpdir(data) as tmp_dir:
            tmp_dir = utils.safe_makedir(os.path.join(tmp_dir, dd.get_sample_name(data)))
            gender_chroms = [x.name for x in ref.file_contigs(dd.get_ref_file(data)) if chromhacks.is_sex(x.name)]
            gender_args = "--sex %s" % (",".join(gender_chroms)) if gender_chroms else ""
            cmd = "goleft indexcov --directory {tmp_dir} {gender_args} -- {bam_file}"
            try:
                do.run(cmd.format(**locals()), "QC: goleft indexcov")
            except subprocess.CalledProcessError as msg:
                if not ("indexcov: no usable" in str(msg) or
                        ("indexcov: expected" in str(msg) and "sex chromosomes, found:" in str(msg))):
                    raise
            for out_file in out_files:
                orig_file = os.path.join(tmp_dir, os.path.basename(out_file))
                if utils.file_exists(orig_file):
                    utils.copy_plus(orig_file, out_file)
    # MultiQC needs non-gzipped/BED inputs so unpack the file
    out_bed = out_files[-1].replace(".bed.gz", ".tsv")
    if utils.file_exists(out_files[-1]) and not utils.file_exists(out_bed):
        with transaction.file_transaction(data, out_bed) as tx_out_bed:
            cmd = "gunzip -c %s > %s" % (out_files[-1], tx_out_bed)
            do.run(cmd, "Unpack indexcov BED file")
    out_files[-1] = out_bed
    return [x for x in out_files if utils.file_exists(x)]
