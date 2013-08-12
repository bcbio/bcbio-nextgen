"""Provide piped, no disk-IO, BAM preparation for variant calling.
Handles independent analysis of chromosome regions, allowing parallel
runs of this step.
"""
import os

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils, shared
from bcbio.provenance import do
from bcbio.variation import realign, recalibrate

# ## GATK/Picard preparation

def region_to_gatk(region):
    if isinstance(region, (list, tuple)):
        chrom, start, end = region
        return "%s:%s-%s" % (chrom, start + 1, end)
    else:
        return region

def _gatk_extract_reads_cl(data, region, prep_params, tmp_dir):
    """Use GATK to extract reads from full BAM file, recalibrating if configured.
    """
    broad_runner = broad.runner_from_config(data["config"])
    args = ["-T", "PrintReads",
            "-L", region_to_gatk(region),
            "-R", data["sam_ref"],
            "-I", data["work_bam"]]
    if prep_params["recal"] == "gatk":
        if _recal_has_reads(data["prep_recal"]):
            args += ["-BQSR", data["prep_recal"]]
    elif prep_params["recal"]:
        raise NotImplementedError("Recalibration method %s" % prep_params["recal"])
    return broad_runner.cl_gatk(args, tmp_dir)

def _recal_has_reads(in_file):
    with open(in_file) as in_handle:
        return not in_handle.readline().startswith("# No aligned reads")

def _piped_input_cl(data, region, tmp_dir, out_base_file, prep_params):
    """Retrieve the commandline for streaming input into preparation step.
    If marking duplicates, this requires writing an intermediate file since
    MarkDuplicates uses multiple passed on an input.
    """
    broad_runner = broad.runner_from_config(data["config"])
    cl = _gatk_extract_reads_cl(data, region, prep_params, tmp_dir)
    if prep_params["dup"] == "picard":
        sel_file = "%s-select%s" % os.path.splitext(out_base_file)
        if not utils.file_exists(sel_file):
            with file_transaction(sel_file) as tx_out_file:
                cl += ["-o", tx_out_file]
                do.run(cl, "GATK: PrintReads {0}".format(region), data)
        dup_metrics = "%s-dup.dup_metrics" % os.path.splitext(out_base_file)[0]
        compression = "5" if prep_params["realign"] == "gatk" else "0"
        cl = broad_runner.cl_picard("MarkDuplicates",
                                    [("INPUT", sel_file),
                                     ("OUTPUT", "/dev/stdout"),
                                     ("METRICS_FILE", dup_metrics),
                                     ("PROGRAM_RECORD_ID", "null"),
                                     ("COMPRESSION_LEVEL", compression),
                                     ("TMP_DIR", tmp_dir)])
    elif not prep_params["dup"]:
        sel_file = data["work_bam"]
    else:
        raise ValueError("Duplication approach not supported with GATK: %s" % prep_params["dup"])
    broad_runner.run_fn("picard_index", sel_file)
    return sel_file, " ".join(cl)

def _piped_realign_gatk(data, region, cl, out_base_file, tmp_dir, prep_params):
    """Perform realignment with GATK, using input commandline.
    GATK requires writing to disk and indexing before realignment.
    """
    broad_runner = broad.runner_from_config(data["config"])
    pa_bam = "%s-prealign%s" % os.path.splitext(out_base_file)
    if not utils.file_exists(pa_bam):
        with file_transaction(pa_bam) as tx_out_file:
            pipe = ">" if prep_params["dup"] else "-o"
            cmd = "{cl} {pipe} {tx_out_file}".format(**locals())
            do.run(cmd, "GATK pre-alignment {0}".format(region), data)
    broad_runner.run_fn("picard_index", pa_bam)
    dbsnp_vcf = shared.configured_ref_file("dbsnp", data["config"], data["sam_ref"])
    recal_file = realign.gatk_realigner_targets(broad_runner, pa_bam, data["sam_ref"],
                                                dbsnp=dbsnp_vcf, region=region_to_gatk(region))
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

def _piped_bamprep_region_gatk(data, region, prep_params, out_file, tmp_dir):
    """Perform semi-piped BAM preparation using Picard/GATK tools.
    """
    broad_runner = broad.runner_from_config(data["config"])
    cur_bam, cl = _piped_input_cl(data, region, tmp_dir, out_file, prep_params)
    if not prep_params["realign"]:
        prerecal_bam = None
    elif prep_params["realign"] == "gatk":
        prerecal_bam, cl = _piped_realign_gatk(data, region, cl, out_file, tmp_dir,
                                               prep_params)
    else:
        raise NotImplementedError("Realignment method: %s" % prep_params["realign"])
    with file_transaction(out_file) as tx_out_file:
        out_flag = ("-o" if prep_params["realign"] == "gatk"
                    or (not prep_params["realign"] and not prep_params["dup"])
                    else ">")
        cmd = "{cl} {out_flag} {tx_out_file}".format(**locals())
        do.run(cmd, "GATK: realign {0}".format(region), data)
        _cleanup_tempfiles(data, [cur_bam, prerecal_bam])

# ## Full-piped approaches

def _piped_dedup_recal_cmd(data, prep_params, tmp_dir):
    """Generate de-duplication and recalibration commandline.
    """
    if prep_params["dup"] == "bamutil":
        assert prep_params["recal"] in ["bamutil", False], \
            "Cannot handle recalibration approach %s with bamutil dedup" % prep_params["recal"]
        out_stream = "-.ubam" if prep_params["realign"] else "-.bam"
        return "| " + recalibrate.bamutil_dedup_recal_cl("-.ubam", out_stream, data,
                                                         prep_params["recal"] == "bamutil")
    elif prep_params["dup"] == "samtools":
        samtools = config_utils.get_program("samtools", data["config"])
        return "| " + "{samtools} rmdup - -".format(**locals())
    elif prep_params["dup"]:
        raise ValueError("Unexpected deduplication approach: %s" % prep_params["dup"])
    else:
        return ""

def _piped_realign_cmd(data, prep_params, tmp_dir):
    """Generate piped realignment commandline.
    """
    if prep_params["realign"] == "gkno":
        return "| " + realign.gkno_realigner_cl(data["sam_ref"], data["config"])
    elif prep_params["realign"]:
        raise ValueError("Unexpected realignment approach: %s" % prep_params["realign"])
    else:
        return ""

def _piped_extract_recal_cmd(data, region, prep_params, tmp_dir):
    """Extract region to process, potentially performing recalibration.
    Combines extraction and recalibration supported by GATK.
    """
    config = data["config"]
    samtools = config_utils.get_program("samtools", config)
    out_type = "-u" if prep_params["dup"] or prep_params["realign"] else "-b"
    if not prep_params.get("recal"):
        prep_region = region_to_gatk(region)
        in_file = data["work_bam"]
        cmd = "{samtools} view {out_type} {in_file} {prep_region}"
        return cmd.format(**locals())
    elif prep_params["recal"] == "gatk":
        cl = _gatk_extract_reads_cl(data, region, prep_params, tmp_dir)
        cl += ["--logging_level", "ERROR"]
        cmd = "{samtools} view -S {out_type} -"
        return " ".join(cl) + " | " + cmd.format(**locals())
    else:
        raise ValueError("Unexpected recalibration approach: %s" % prep_params["recal"])

def _piped_bamprep_region_fullpipe(data, region, prep_params, out_file, tmp_dir):
    """Perform fully piped BAM preparation using non-GATK/Picard tools.
    """
    with file_transaction(out_file) as tx_out_file:
        extract_recal_cmd = _piped_extract_recal_cmd(data, region, prep_params, tmp_dir)
        dedup_cmd = _piped_dedup_recal_cmd(data, prep_params, tmp_dir)
        realign_cmd = _piped_realign_cmd(data, prep_params, tmp_dir)
        cmd = "{extract_recal_cmd} {dedup_cmd} {realign_cmd}  > {tx_out_file}"
        cmd = cmd.format(**locals())
        do.run(cmd, "Piped post-alignment bamprep {0}".format(region), data)

# ## Shared functionality

def _get_prep_params(data):
    """Retrieve configuration parameters with defaults for preparing BAM files.
    """
    algorithm = data["config"]["algorithm"]
    dup_param = algorithm.get("mark_duplicates", True)
    dup_param = "picard" if dup_param is True else dup_param
    recal_param = algorithm.get("recalibrate", True)
    recal_param = "gatk" if recal_param is True else recal_param
    realign_param = algorithm.get("realign", True)
    realign_param = "gatk" if realign_param is True else realign_param
    all_params = set([dup_param, realign_param])
    return {"dup": dup_param, "recal": recal_param, "realign": realign_param,
            "all_pipe": "gatk" not in all_params and "picard" not in all_params}

def _piped_bamprep_region(data, region, out_file, tmp_dir):
    """Do work of preparing BAM input file on the selected region.
    """
    prep_params = _get_prep_params(data)
    if prep_params["all_pipe"]:
        _piped_bamprep_region_fullpipe(data, region, prep_params, out_file, tmp_dir)
    else:
        _piped_bamprep_region_gatk(data, region, prep_params, out_file, tmp_dir)

def piped_bamprep(data, region=None, out_file=None):
    """Perform full BAM preparation using pipes to avoid intermediate disk IO.

    Handles de-duplication, recalibration and realignment of original BAMs.
    """
    utils.safe_makedir(os.path.dirname(out_file))
    if region[0] == "nochrom":
        prep_bam = shared.write_nochr_reads(data["work_bam"], out_file)
    elif region[0] == "noanalysis":
        prep_bam = shared.write_noanalysis_reads(data["work_bam"], region[1], out_file,
                                                 data["config"])
    else:
        if not utils.file_exists(out_file):
            with utils.curdir_tmpdir() as tmp_dir:
                _piped_bamprep_region(data, region, out_file, tmp_dir)
        prep_bam = out_file
    broad_runner = broad.runner_from_config(data["config"])
    broad_runner.run_fn("picard_index", prep_bam)
    data["work_bam"] = prep_bam
    data["region"] = region
    return [data]
