"""Provide piped, no disk-IO, BAM preparation for variant calling.
Handles independent analysis of chromosome regions, allowing parallel
runs of this step.
"""
import os

from bcbio import bam, broad, utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.pipeline import shared
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import realign

# ## GATK/Picard preparation

def region_to_gatk(region):
    if isinstance(region, (list, tuple)):
        chrom, start, end = region
        return "%s:%s-%s" % (chrom, start + 1, end)
    else:
        return region

def _gatk_extract_reads_cl(data, region, prep_params, tmp_dir):
    """Use GATK to extract reads from full BAM file.
    """
    args = ["PrintReads",
            "-L", region_to_gatk(region),
            "-R", dd.get_ref_file(data),
            "-I", data["work_bam"]]
    # GATK3 back compatibility, need to specify analysis type
    if "gatk4" in dd.get_tools_off(data):
        args = ["--analysis_type"] + args
    runner = broad.runner_from_config(data["config"])
    return runner.cl_gatk(args, tmp_dir)

def _piped_input_cl(data, region, tmp_dir, out_base_file, prep_params):
    """Retrieve the commandline for streaming input into preparation step.
    """
    return data["work_bam"], _gatk_extract_reads_cl(data, region, prep_params, tmp_dir)

def _piped_realign_gatk(data, region, cl, out_base_file, tmp_dir, prep_params):
    """Perform realignment with GATK, using input commandline.
    GATK requires writing to disk and indexing before realignment.
    """
    broad_runner = broad.runner_from_config(data["config"])
    pa_bam = "%s-prealign%s" % os.path.splitext(out_base_file)
    if not utils.file_exists(pa_bam):
        with file_transaction(data, pa_bam) as tx_out_file:
            cmd = "{cl} -o {tx_out_file}".format(**locals())
            do.run(cmd, "GATK re-alignment {0}".format(region), data)
    bam.index(pa_bam, data["config"])
    realn_file = realign.gatk_realigner_targets(broad_runner, pa_bam, dd.get_ref_file(data), data["config"],
                                                region=region_to_gatk(region),
                                                known_vrns=dd.get_variation_resources(data))
    realn_cl = realign.gatk_indel_realignment_cl(broad_runner, pa_bam, dd.get_ref_file(data),
                                                 realn_file, tmp_dir, region=region_to_gatk(region),
                                                 known_vrns=dd.get_variation_resources(data))
    return pa_bam, realn_cl

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
    with file_transaction(data, out_file) as tx_out_file:
        out_flag = ("-o" if (prep_params["realign"] == "gatk"
                             or not prep_params["realign"])
                    else ">")
        cmd = "{cl} {out_flag} {tx_out_file}".format(**locals())
        do.run(cmd, "GATK: realign {0}".format(region), data)
        _cleanup_tempfiles(data, [cur_bam, prerecal_bam])

# ## Shared functionality

def _get_prep_params(data):
    """Retrieve configuration parameters with defaults for preparing BAM files.
    """
    realign_param = dd.get_realign(data)
    realign_param = "gatk" if realign_param is True else realign_param
    return {"realign": realign_param}

def _need_prep(data):
    prep_params = _get_prep_params(data)
    return prep_params["realign"]

def _piped_bamprep_region(data, region, out_file, tmp_dir):
    """Do work of preparing BAM input file on the selected region.
    """
    if _need_prep(data):
        prep_params = _get_prep_params(data)
        _piped_bamprep_region_gatk(data, region, prep_params, out_file, tmp_dir)
    else:
        raise ValueError("No realignment specified")

def piped_bamprep(data, region=None, out_file=None):
    """Perform full BAM preparation using pipes to avoid intermediate disk IO.

    Handles realignment of original BAMs.
    """
    data["region"] = region
    if not _need_prep(data):
        return [data]
    else:
        utils.safe_makedir(os.path.dirname(out_file))
        if region[0] == "nochrom":
            prep_bam = shared.write_nochr_reads(data["work_bam"], out_file, data["config"])
        elif region[0] == "noanalysis":
            prep_bam = shared.write_noanalysis_reads(data["work_bam"], region[1], out_file,
                                                     data["config"])
        else:
            if not utils.file_exists(out_file):
                with tx_tmpdir(data) as tmp_dir:
                    _piped_bamprep_region(data, region, out_file, tmp_dir)
            prep_bam = out_file
        bam.index(prep_bam, data["config"])
        data["work_bam"] = prep_bam
        return [data]
