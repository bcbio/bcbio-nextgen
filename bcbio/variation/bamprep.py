"""Provide piped, no disk-IO, BAM preparation for variant calling.
Handles independent analysis of chromosome regions, allowing parallel
runs of this step.
"""
import os
import subprocess

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import shared
from bcbio.variation import realign

def region_to_gatk(region):
    if isinstance(region, (list, tuple)):
        chrom, start, end = region
        return "%s:%s-%s" % (chrom, start + 1, end)
    else:
        return region

def _gatk_extract_reads_cl(data, region, tmp_dir):
    """Use GATK to extract reads from full BAM file, recalibrating if configured.
    """
    broad_runner = broad.runner_from_config(data["config"])
    algorithm = data["config"]["algorithm"]
    args = ["-T", "PrintReads",
            "-L", region_to_gatk(region),
            "-R", data["sam_ref"],
            "-I", data["work_bam"]]
    recal_config = algorithm.get("recalibrate", True)
    if recal_config in ["gatk", True]:
        args += ["-BQSR", data["prep_recal"]]
    elif recal_config:
        raise NotImplementedError("Recalibration method %s" %  recal_config)
    return broad_runner.cl_gatk(args, tmp_dir)

def _piped_input_cl(data, region, tmp_dir, out_base_file, realign_param):
    """Retrieve the commandline for streaming input into preparation step.
    If marking duplicates, this requires writing an intermediate file since
    MarkDuplicates uses multiple passed on an input.
    """
    broad_runner = broad.runner_from_config(data["config"])
    cl = _gatk_extract_reads_cl(data, region, tmp_dir)
    algorithm = data["config"]["algorithm"]
    if algorithm.get("mark_duplicates", True):
        sel_file = "%s-select%s" % os.path.splitext(out_base_file)
        if not utils.file_exists(sel_file):
            with file_transaction(sel_file) as tx_out_file:
                cl += ["-o", tx_out_file]
                subprocess.check_call(cl)
        dup_metrics = "%s-dup.dup_metrics" % os.path.splitext(out_base_file)[0]
        compression = "5" if realign_param == "gatk" else "0"
        cl = broad_runner.cl_picard("MarkDuplicates",
                                    [("INPUT", sel_file),
                                     ("OUTPUT", "/dev/stdout"),
                                     ("METRICS_FILE", dup_metrics),
                                     ("PROGRAM_RECORD_ID", "null"),
                                     ("COMPRESSION_LEVEL", compression),
                                     ("TMP_DIR", tmp_dir)])
    else:
        sel_file = data["work_bam"]
    broad_runner.run_fn("picard_index", sel_file)
    return sel_file, " ".join(cl)

def _piped_realign_gatk(data, region, cl, out_base_file, tmp_dir):
    """Perform realignment with GATK, using input commandline.
    GATK requires writing to disk and indexing before realignment.
    """
    broad_runner = broad.runner_from_config(data["config"])
    pa_bam = "%s-prealign%s" % os.path.splitext(out_base_file)
    if not utils.file_exists(pa_bam):
        with file_transaction(pa_bam) as tx_out_file:
            subprocess.check_call("{cl} > {tx_out_file}".format(**locals()), shell=True)
    broad_runner.run_fn("picard_index", pa_bam)
    recal_file = realign.gatk_realigner_targets(broad_runner, pa_bam, data["sam_ref"],
                      dbsnp=shared.configured_ref_file("dbsnp", data["config"], data["sam_ref"]),
                      region=region_to_gatk(region))
    recal_cl = realign.gatk_indel_realignment_cl(broad_runner, pa_bam, data["sam_ref"],
                                                 recal_file, tmp_dir, region=region_to_gatk(region))
    return pa_bam, " ".join(recal_cl)

def _cleanup_tempfiles(data, tmp_files):
    for tmp_file in tmp_files:
        if tmp_file and tmp_file != data["work_bam"]:
            for ext in [".bam", ".bam.bai", ".bai"]:
                fname = "%s%s" % (os.path.splitext(tmp_file)[0], ext)
                if os.path.exists(fname):
                    os.remove(fname)

def _piped_bamprep_region(data, region, out_file, tmp_dir):
    """Do work of preparing BAM input file on the selected region.
    """
    broad_runner = broad.runner_from_config(data["config"])
    algorithm = data["config"]["algorithm"]
    realign_param = algorithm.get("realign", "gatk")
    realign_param = "gatk" if realign_param is True else realign_param
    cur_bam, cl = _piped_input_cl(data, region, tmp_dir, out_file, realign_param)
    if not realign_param:
        prerecal_bam = None
    elif realign_param == "gatk":
        prerecal_bam, cl = _piped_realign_gatk(data, region, cl, out_file, tmp_dir)
    else:
        raise NotImplementedError("Realignment method: %s" % realign_param)
    with file_transaction(out_file) as tx_out_file:
        out_flag = "-o" if realign_param == "gatk" else ">"
        subprocess.check_call("{cl} {out_flag} {tx_out_file}".format(**locals()), shell=True)
        _cleanup_tempfiles(data, [cur_bam, prerecal_bam])

def piped_bamprep(data, region=None, out_file=None):
    """Perform full BAM preparation using pipes to avoid intermediate disk IO.

    Handles de-duplication, recalibration and realignment of original BAMs.
    """
    utils.safe_makedir(os.path.dirname(out_file))
    if region[0] == "nochrom":
        prep_bam = shared.write_nochr_reads(data["work_bam"], out_file)
    elif region[0] == "noanalysis":
        prep_bam = shared.write_noanalysis_reads(data["work_bam"], region[1], out_file)
    else:
        if not utils.file_exists(out_file):
            with utils.curdir_tmpdir() as tmp_dir:
                _piped_bamprep_region(data, region, out_file, tmp_dir)
        prep_bam = out_file
    data["work_bam"] = prep_bam
    data["region"] = region
    return [data]
