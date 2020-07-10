import os
import subprocess
import glob

from bcbio import utils
from bcbio.provenance import do
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio import bam
from bcbio.chipseq import antibodies
from bcbio.log import logger
from bcbio.distributed.transaction import file_transaction

def run(name, chip_bam, input_bam, genome_build, out_dir, method, resources, data):
    """
    Run macs2 for chip and input samples avoiding
    errors due to samples.
    """
    # output file name need to have the caller name
    config = dd.get_config(data)
    out_file = os.path.join(out_dir, name + "_peaks_macs2.xls")
    macs2_file = os.path.join(out_dir, name + "_peaks.xls")
    if utils.file_exists(out_file):
        _compress_and_sort_bdg_files(out_dir, data)
        return _get_output_files(out_dir)
    macs2 = config_utils.get_program("macs2", config)
    antibody = antibodies.ANTIBODIES.get(dd.get_antibody(data).lower(), None)
    if antibody:
        logger.info(f"{antibody.name} specified, using {antibody.peaktype} peak settings.")
        peaksettings = select_peak_parameters(antibody)
    elif method == "atac":
        logger.info(f"ATAC-seq specified, using narrow peak settings.")
        peaksettings = " "
    else:
        peaksettings = " "
    options = " ".join(resources.get("macs2", {}).get("options", ""))
    genome_size = bam.fasta.total_sequence_length(dd.get_ref_file(data))
    genome_size = "" if options.find("-g") > -1 else "-g %s" % genome_size
    paired = "-f BAMPE" if bam.is_paired(chip_bam) else ""
    with utils.chdir(out_dir):
        cmd = _macs2_cmd(data)
        cmd += peaksettings
        try:
            do.run(cmd.format(**locals()), "macs2 for %s" % name)
            utils.move_safe(macs2_file, out_file)
        except subprocess.CalledProcessError:
            raise RuntimeWarning("macs2 terminated with an error. "
                                 "Please, check the message and report "
                                 "error if it is related to bcbio. "
                                 "You can add specific options for the sample "
                                 "setting resources as explained in docs: "
                                 "https://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#sample-specific-resources")
    _compress_and_sort_bdg_files(out_dir, data)
    return _get_output_files(out_dir)

def _get_output_files(out_dir):
    fns = [os.path.abspath(fn) for fn in glob.glob(os.path.join(out_dir, "*"))]
    peaks = None
    for fn in fns:
        if fn.endswith("narrowPeak"):
            peaks = fn
            break
        elif fn.endswith("broadPeak"):
            peaks = fn
            break
    return {"main": peaks, "macs2": fns}

def _compress_and_sort_bdg_files(out_dir, data):
    for fn in glob.glob(os.path.join(out_dir, "*bdg")):
        out_file = fn + ".gz"
        if utils.file_exists(out_file):
            continue
        bedtools = config_utils.get_program("bedtools", data)
        with file_transaction(out_file) as tx_out_file:
            cmd = f"sort -k1,1 -k2,2n {fn} | bgzip -c > {tx_out_file}"
            message = f"Compressing and sorting {fn}."
            do.run(cmd, message)

def _macs2_cmd(data):
    """Main command for macs2 tool."""
    method = dd.get_chip_method(data)
    if method.lower() == "chip":
        cmd = ("{macs2} callpeak -t {chip_bam} -c {input_bam} {paired} "
               "{genome_size} -n {name} --bdg {options} ")
    elif method.lower() == "atac":
        cmd = ("{macs2} callpeak -t {chip_bam} --nomodel "
               " {paired} {genome_size} -n {name} --bdg {options}"
               " --nolambda --keep-dup all")
    else:
        raise ValueError("chip_method should be chip or atac.")
    return cmd

def select_peak_parameters(antibody):
    if antibody.peaktype == "broad":
        return " --broad --broad-cutoff 0.05"
    elif antibody.peaktype == "narrow":
        return ""
    else:
        raise ValueError(f"{antibody.peaktype} not recognized.")
