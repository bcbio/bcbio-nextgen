import os
import toolz as tz
from collections import namedtuple

from bcbio import utils
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.log import logger
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio import bam
from bcbio.rnaseq import gtf
from bcbio.heterogeneity import chromhacks

# ranges taken from Buenrostro, Nat. Methods 10, 1213–1218 (2013).
ATACRange = namedtuple('ATACRange', ['label', 'min', 'max'])
ATACRanges = {"NF": ATACRange("NF", 0, 100),
              "MN": ATACRange("MN", 180, 247),
              "DN": ATACRange("DN", 315, 473),
              "TN": ATACRange("TN", 558, 615)}

def calculate_complexity_metrics(work_bam, data):
    """
    the work_bam should have duplicates marked but not removed
    mitochondrial reads should be removed 
    """
    bedtools = config_utils.get_program("bedtools", dd.get_config(data))
    work_dir = dd.get_work_dir(data)
    metrics_dir = os.path.join(work_dir, "metrics", "atac")
    utils.safe_makedir(metrics_dir)
    metrics_file = os.path.join(metrics_dir,
                                f"{dd.get_sample_name(data)}-atac-metrics.csv")
    if utils.file_exists(metrics_file):
        data = tz.assoc_in(data, ['atac', 'complexity_metrics_file'], metrics_file)
        return data
    # BAM file must be sorted by read name
    work_bam = bam.sort(work_bam, dd.get_config(data), order="queryname")
    with file_transaction(metrics_file) as tx_metrics_file:
        with open(tx_metrics_file, "w") as out_handle:
            out_handle.write("mt,m0,m1,m2\n")
        if bam.is_paired(work_bam):
            cmd = (f"{bedtools} bamtobed -bedpe -i {work_bam} | "
                "awk 'BEGIN{OFS=\"\\t\"}{print $1,$2,$4,$6,$9,$10}' | "
                "sort | "
                "uniq -c | "
                "awk 'BEGIN{mt=0;m0=0;m1=0;m2=0}($1==1){m1=m1+1} ($1==2){m2=m2+1}{m0=m0+1}{mt=mt+$1} END{printf \"%d,%d,%d,%d\\n\", mt,m0,m1,m2}' >> "
                f"{tx_metrics_file}")
        else:
            cmd = (f"{bedtools} bamtobed -i {work_bam} | "
                   "awk 'BEGIN{OFS=\"\\t\"}{print $1,$2,$3,$6}' | "
                   "sort | "
                   "uniq -c | "
                   "awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1}{m0=m0+1}{mt=mt+$1} END{printf \"%d,%d,%d,%d\\n\", mt,m0,m1,m2}' >> "
                   f"{tx_metrics_file}")
        message = f"Calculating ATAC-seq complexity metrics on {work_bam}, saving as {metrics_file}."
        do.run(cmd, message)
    data = tz.assoc_in(data, ['atac', 'complexity_metrics_file'], metrics_file)
    return data

def calculate_encode_complexity_metrics(data):
    metrics_file = tz.get_in(['atac', 'complexity_metrics_file'], data, None)
    if not metrics_file:
        return {}
    else:
        with open(metrics_file) as in_handle:
            header = next(in_handle).strip().split(",")
            values = next(in_handle).strip().split(",")
    raw_metrics = {h: int(v) for h, v in zip(header, values)}
    metrics = {"PBC1": raw_metrics["m1"] / raw_metrics["m0"],
               "NRF": raw_metrics["m0"] / raw_metrics["mt"]}
    if raw_metrics["m2"] == 0:
        PBC2 = 0
    else:
        PBC2 = raw_metrics["m1"] / raw_metrics["m2"]
    metrics["PBC2"] = PBC2
    metrics["bottlenecking"] = get_bottlenecking_flag(metrics["PBC1"], metrics["PBC2"])
    metrics["complexity"] = get_complexity_flag(metrics["NRF"])
    return(metrics)

def get_bottlenecking_flag(PBC1, PBC2):
    if PBC1 < 0.7 or PBC2 < 1:
        return "severe"
    elif PBC1 <= 0.9 or PBC2 <= 3:
        return "moderate"
    else:
        return "none"

def get_complexity_flag(NRF):
    if NRF < 0.7:
        return "concerning"
    elif NRF < 0.9:
        return "acceptable"
    else:
        return "ideal"


def split_ATAC(data, bam_file=None):
    """
    splits a BAM into nucleosome-free (NF) and mono/di/tri nucleosome BAMs based
    on the estimated insert sizes
    uses the current working BAM file if no BAM file is supplied
    """
    sambamba = config_utils.get_program("sambamba", data)
    num_cores = dd.get_num_cores(data)
    base_cmd = f'{sambamba} view --format bam --nthreads {num_cores} '
    bam_file = bam_file if bam_file else dd.get_work_bam(data)
    out_stem = os.path.splitext(bam_file)[0]
    split_files = {}
    # we can only split these fractions from paired runs
    if not bam.is_paired(bam_file):
        split_files["full"] = bam_file
        data = tz.assoc_in(data, ['atac', 'align'], split_files)
        return data
    for arange in ATACRanges.values():
        out_file = f"{out_stem}-{arange.label}.bam"
        if not utils.file_exists(out_file):
            with file_transaction(out_file) as tx_out_file:
                cmd = base_cmd +\
                    f'-F "template_length > {arange.min} and template_length < {arange.max}" ' +\
                    f'{bam_file} > {tx_out_file}'
                message = f'Splitting {arange.label} regions from {bam_file}.'
                do.run(cmd, message)
            bam.index(out_file, dd.get_config(data))
        split_files[arange.label] = out_file
    split_files["full"] = bam_file
    data = tz.assoc_in(data, ['atac', 'align'], split_files)
    return data

def run_ataqv(data):
    if not dd.get_chip_method(data) == "atac":
        return None
    work_dir = dd.get_work_dir(data)
    sample_name = dd.get_sample_name(data)
    out_dir = os.path.join(work_dir, "qc", sample_name, "ataqv")
    peak_file = get_full_peaks(data)
    bam_file = get_unfiltered_bam(data)
    out_file = os.path.join(out_dir, sample_name + ".ataqv.json.gz")
    if not peak_file:
        logger.info(f"Full peak file for {sample_name} not found, skipping ataqv")
        return None
    if not bam_file:
        logger.info(f"Unfiltered BAM file for {sample_name} not found, skipping ataqv")
        return None
    if utils.file_exists(out_file):
        return out_file
    tss_bed_file = os.path.join(out_dir, "TSS.bed")
    tss_bed_file = gtf.get_tss_bed(dd.get_gtf_file(data), tss_bed_file, data, padding=0)
    if chromhacks.is_human(data):
        organism = "human"
        autosomal_reference_flag = ""
    elif chromhacks.is_mouse(data):
        organism = "mouse"
        autosomal_reference_flag = ""
    else:
        autosomal_reference = os.path.join(out_dir, "autosomal.txt")
        autosomal_reference = _make_autosomal_reference_file(autosomal_reference, data)
        organism = "None"
        autosomal_reference_flag = f"--autosomal-reference-file {autosomal_reference} "
    ataqv = config_utils.get_program("ataqv", data)
    mitoname = chromhacks.get_mitochondrial_chroms(data)[0]
    if not ataqv:
        logger.info(f"ataqv executable not found, skipping running ataqv.")
        return None
    with file_transaction(out_file) as tx_out_file:
        cmd = (f"{ataqv} --peak-file {peak_file} --name {sample_name} --metrics-file {tx_out_file} "
               f"--tss-file {tss_bed_file} {autosomal_reference_flag} "
               f"--ignore-read-groups --mitochondrial-reference-name {mitoname} "
               f"--tss-extension 1000 "
               f"{organism} {bam_file}")
        message = f"Running ataqv on {sample_name}."
        do.run(cmd, message)
    return out_file

def _make_autosomal_reference_file(out_file, data):
    """
    for many organisms we don't know in bcbio what chromosomes are what, for now include
    everything non-mitochondrial
    """
    if utils.file_exists(out_file):
        return out_file
    nonmito = chromhacks.get_nonmitochondrial_chroms(data)
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            for chrom in nonmito:
                print(f"{chrom}", file=out_handle)
    return out_file

def get_NF_bam(data):
    """
    get the nucleosome free BAM file for ATAC-seq if it exists
    """
    return tz.get_in(("atac", "align", "NF"), data, None)

def get_NF_peaks(data):
    """
    get the nucleosome free peak file for ATAC-seq if it exists
    """
    peak_files = tz.get_in(("peaks_files", "NF", "macs2"), data, [])
    for f in peak_files:
        if f.endswith("narrowPeak") or f.endswith("broadPeak"):
            return f
    return None

def get_unfiltered_bam(data):
    """
    get the nucleosome free BAM file for ATAC-seq if it exists
    """
    return tz.get_in(("chipseq", "align", "unfiltered"), data, None)

def get_full_peaks(data):
    """
    get the nucleosome free peak file for ATAC-seq if it exists
    """
    peak_files = tz.get_in(("peaks_files", "full", "macs2"), data, [])
    for f in peak_files:
        if f.endswith("narrowPeak") or f.endswith("broadPeak"):
            return f
    return None

def create_ataqv_report(samples):
    """
    make the ataqv report from a set of ATAC-seq samples
    """
    data = samples[0][0]
    new_samples = []
    reportdir = os.path.join(dd.get_work_dir(data), "qc", "ataqv")
    sentinel = os.path.join(reportdir, "index.html")
    if utils.file_exists(sentinel):
        ataqv_output = {"base": sentinel, "secondary": get_ataqv_report_files(reportdir)}
        new_data = []
        for data in dd.sample_data_iterator(samples):
            data = tz.assoc_in(data, ["ataqv_report"], ataqv_output)
            new_data.append(data)
        return dd.get_samples_from_datalist(new_data)
    mkarv = config_utils.get_program("mkarv", dd.get_config(data))
    ataqv_files = []
    for data in dd.sample_data_iterator(samples):
        qc = dd.get_summary_qc(data)
        ataqv_file = tz.get_in(("ataqv", "base"), qc, None)
        if ataqv_file and utils.file_exists(ataqv_file):
            ataqv_files.append(ataqv_file)
    if not ataqv_files:
        return samples
    ataqv_json_file_string = " ".join(ataqv_files)
    with file_transaction(reportdir) as txreportdir:
        cmd = f"{mkarv} {txreportdir} {ataqv_json_file_string}"
        message = f"Creating ataqv report from {ataqv_json_file_string}."
        do.run(cmd, message)
    new_data = []
    ataqv_output = {"base": sentinel, "secondary": get_ataqv_report_files(reportdir)}
    for data in dd.sample_data_iterator(samples):
        data = tz.assoc_in(data, ["ataqv_report"], ataqv_output)
        new_data.append(data)
    return dd.get_samples_from_datalist(new_data)

def get_ataqv_report_files(reportdir):
    files = []
    for r, d, f in os.walk(reportdir):
        for file in f:
            f = os.path.join(r, file)
            if utils.file_exists(f):
                files.append(f)
    return files
