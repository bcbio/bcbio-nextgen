"""Coverage based QC calculations.
"""
import os
import subprocess

from bcbio.bam import ref, sambamba, utils
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

    if dd.get_coverage(data) and dd.get_coverage(data) not in ["None"]:
        cov_bed_file = bedutils.clean_file(dd.get_coverage(data), data, prefix="cov-", simple=True)
        merged_bed_file = bedutils.merge_overlaps(cov_bed_file, data)
        target_name = "coverage"
    elif dd.get_coverage_interval(data) != "genome":
        merged_bed_file = dd.get_variant_regions_merged(data)
        target_name = "variant_regions"
    else:
        merged_bed_file = None
        target_name = "genome"

    avg_depth = cov.get_average_coverage(data, bam_file, merged_bed_file, target_name)
    out['Avg_coverage'] = avg_depth

    samtools_stats_dir = os.path.join(out_dir, os.path.pardir, out_dir)
    from bcbio.qc import samtools
    samtools_stats = samtools.run(bam_file, data, samtools_stats_dir)

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
        mapped_unique = sambamba.number_of_mapped_reads(data, bam_file, keep_dups=False)
    out['Mapped_unique_reads'] = mapped_unique

    if merged_bed_file:
        ontarget = sambamba.number_of_mapped_reads(
            data, bam_file, keep_dups=False, bed_file=merged_bed_file, target_name=target_name)
        out["Ontarget_unique_reads"] = ontarget
        if mapped_unique:
            out["Ontarget_pct"] = 100.0 * ontarget / mapped_unique
            out['Offtarget_pct'] = 100.0 * (mapped_unique - ontarget) / mapped_unique
            if dd.get_coverage_interval(data) != "genome":
                # Skip padded calculation for WGS even if the "coverage" file is specified
                # the padded statistic makes only sense for exomes and panels
                padded_bed_file = bedutils.get_padded_bed_file(merged_bed_file, 200, data)
                ontarget_padded = sambamba.number_of_mapped_reads(
                    data, bam_file, keep_dups=False, bed_file=padded_bed_file, target_name=target_name + "_padded")
                out["Ontarget_padded_pct"] = 100.0 * ontarget_padded / mapped_unique
        if total_reads:
            out['Usable_pct'] = 100.0 * ontarget / total_reads

    region_coverage_file = cov.coverage_region_detailed_stats(data, out_dir,
                                                              extra_cutoffs=set([max(1, int(avg_depth * 0.8))]))
    indexcov_files = _goleft_indexcov(bam_file, data, out_dir)
    out_files = [x for x in [region_coverage_file] + indexcov_files if x and utils.file_exists(x)]
    out = {"metrics": out}
    if len(out_files) > 0:
        out["base"] = out_files[0]
        out["secondary"] = out_files[1:]
    return out

def _goleft_indexcov(bam_file, data, out_dir):
    """Use goleft indexcov to estimate coverage distributions using BAM index.
    """
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
                if "indexcov: no usable" not in str(msg):
                    raise
            for out_file in out_files:
                orig_file = os.path.join(tmp_dir, os.path.basename(out_file))
                if utils.file_exists(orig_file):
                    utils.copy_plus(orig_file, out_file)
    return [x for x in out_files if utils.file_exists(out_file)]
