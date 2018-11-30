"""Qsignature: detection of sample mixups.

https://sourceforge.net/p/adamajava/wiki/qSignature/
"""
import os
import shutil
import subprocess
import xml.etree.ElementTree as ET

import pysam
import toolz as tz

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.provenance import do
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils

def run(bam_file, data, out_dir):
    """ Run SignatureGenerator to create normalize vcf that later will be input of qsignature_summary

    :param bam_file: (str) path of the bam_file
    :param data: (list) list containing the all the dictionary
                     for this sample
    :param out_dir: (str) path of the output

    :returns: (string) output normalized vcf file
    """
    qsig = config_utils.get_program("qsignature", data["config"])
    res_qsig = config_utils.get_resources("qsignature", data["config"])
    jvm_opts = " ".join(res_qsig.get("jvm_opts", ["-Xms750m", "-Xmx8g"]))
    if not qsig:
        logger.info("There is no qsignature tool. Skipping...")
        return None

    position = dd.get_qsig_file(data)
    mixup_check = dd.get_mixup_check(data)
    if mixup_check and mixup_check.startswith("qsignature"):
        utils.safe_makedir(out_dir)
        if not position:
            logger.info("There is no qsignature for this species: %s"
                        % tz.get_in(['genome_build'], data))
            return None
        if mixup_check == "qsignature_full":
            down_bam = bam_file
        else:
            down_bam = _slice_bam_chr21(bam_file, data)
            position = _slice_vcf_chr21(position, out_dir)

        out_name = os.path.basename(down_bam).replace("bam", "qsig.vcf")
        out_file = os.path.join(out_dir, out_name)
        log_file = os.path.join(out_dir, "qsig.log")
        cores = dd.get_cores(data)
        base_cmd = ("{qsig} {jvm_opts} "
                    "org.qcmg.sig.SignatureGenerator "
                    "--noOfThreads {cores} "
                    "-log {log_file} -i {position} "
                    "-i {down_bam} ")
        if not os.path.exists(out_file):
            file_qsign_out = "{0}.qsig.vcf".format(down_bam)
            do.run(base_cmd.format(**locals()), "qsignature vcf generation: %s" % dd.get_sample_name(data))
            if os.path.exists(file_qsign_out):
                with file_transaction(data, out_file) as file_txt_out:
                    shutil.move(file_qsign_out, file_txt_out)
            else:
                raise IOError("File doesn't exist %s" % file_qsign_out)
        return out_file
    return None

def summary(*samples):
    """Run SignatureCompareRelatedSimple module from qsignature tool.

    Creates a matrix of pairwise comparison among samples. The
    function will not run if the output exists

    :param samples: list with only one element containing all samples information
    :returns: (dict) with the path of the output to be joined to summary
    """
    warnings, similar = [], []
    qsig = config_utils.get_program("qsignature", samples[0][0]["config"])
    if not qsig:
        return [[]]
    res_qsig = config_utils.get_resources("qsignature", samples[0][0]["config"])
    jvm_opts = " ".join(res_qsig.get("jvm_opts", ["-Xms750m", "-Xmx8g"]))
    work_dir = samples[0][0]["dirs"]["work"]
    count = 0
    for data in samples:
        data = data[0]
        vcf = tz.get_in(["summary", "qc", "qsignature", "base"], data)
        if vcf:
            count += 1
            vcf_name = dd.get_sample_name(data) + ".qsig.vcf"
            out_dir = utils.safe_makedir(os.path.join(work_dir, "qsignature"))
            if not os.path.lexists(os.path.join(out_dir, vcf_name)):
                os.symlink(vcf, os.path.join(out_dir, vcf_name))
    if count > 0:
        qc_out_dir = utils.safe_makedir(os.path.join(work_dir, "qc", "qsignature"))
        out_file = os.path.join(qc_out_dir, "qsignature.xml")
        out_ma_file = os.path.join(qc_out_dir, "qsignature.ma")
        out_warn_file = os.path.join(qc_out_dir, "qsignature.warnings")
        log = os.path.join(work_dir, "qsignature", "qsig-summary.log")
        if not os.path.exists(out_file):
            with file_transaction(samples[0][0], out_file) as file_txt_out:
                base_cmd = ("{qsig} {jvm_opts} "
                            "org.qcmg.sig.SignatureCompareRelatedSimple "
                            "-log {log} -dir {out_dir} "
                            "-o {file_txt_out} ")
                do.run(base_cmd.format(**locals()), "qsignature score calculation")
        error, warnings, similar = _parse_qsignature_output(out_file, out_ma_file,
                                                            out_warn_file, samples[0][0])
        return [{'total samples': count,
                 'similar samples pairs': len(similar),
                 'warnings samples pairs': len(warnings),
                 'error samples': list(error),
                 'out_dir': qc_out_dir}]
    else:
        return []

def get_qsig_multiqc_files(*samples):
    work_dir = samples[0][0]["dirs"]["work"]
    qc_out_dir = utils.safe_makedir(os.path.join(work_dir, "qc", "qsignature"))
    return [os.path.join(qc_out_dir, "qsignature.ma")]

def _parse_qsignature_output(in_file, out_file, warning_file, data):
    """ Parse xml file produced by qsignature

    :param in_file: (str) with the path to the xml file
    :param out_file: (str) with the path to output file
    :param warning_file: (str) with the path to warning file

    :returns: (list) with samples that could be duplicated

    """
    name = {}
    error, warnings, similar = set(), set(), set()
    same, replicate, related = 0, 0.1, 0.18
    mixup_check = dd.get_mixup_check(data)
    if mixup_check == "qsignature_full":
        same, replicate, related = 0, 0.01, 0.061
    with open(in_file, 'r') as in_handle:
        with file_transaction(data, out_file) as out_tx_file:
            with file_transaction(data, warning_file) as warn_tx_file:
                with open(out_tx_file, 'w') as out_handle:
                    with open(warn_tx_file, 'w') as warn_handle:
                        et = ET.parse(in_handle)
                        for i in list(et.iter('file')):
                            name[i.attrib['id']] = os.path.basename(i.attrib['name']).replace(".qsig.vcf", "")
                        for i in list(et.iter('comparison')):
                            msg = None
                            pair = "-".join([name[i.attrib['file1']], name[i.attrib['file2']]])
                            out_handle.write("%s\t%s\t%s\n" %
                                             (name[i.attrib['file1']], name[i.attrib['file2']], i.attrib['score']))
                            if float(i.attrib['score']) == same:
                                msg = 'qsignature ERROR: read same samples:%s\n'
                                error.add(pair)
                            elif float(i.attrib['score']) < replicate:
                                msg = 'qsignature WARNING: read similar/replicate samples:%s\n'
                                warnings.add(pair)
                            elif float(i.attrib['score']) < related:
                                msg = 'qsignature NOTE: read relative samples:%s\n'
                                similar.add(pair)
                            if msg:
                                logger.info(msg % pair)
                                warn_handle.write(msg % pair)
    return error, warnings, similar

def _slice_bam_chr21(in_bam, data):
    """
    return only one BAM file with only chromosome 21
    """
    sambamba = config_utils.get_program("sambamba", data["config"])
    out_file = "%s-chr%s" % os.path.splitext(in_bam)
    if not utils.file_exists(out_file):
        bam.index(in_bam, data['config'])
        with pysam.Samfile(in_bam, "rb") as bamfile:
            bam_contigs = [c["SN"] for c in bamfile.header["SQ"]]
        chromosome = "21"
        if "chr21" in bam_contigs:
            chromosome = "chr21"
        with file_transaction(data, out_file) as tx_out_file:
            cmd = ("{sambamba} slice -o {tx_out_file} {in_bam} {chromosome}").format(**locals())
            out = subprocess.check_output(cmd, shell=True)
    return out_file

def _slice_vcf_chr21(vcf_file, out_dir):
    """
    Slice chr21 of qsignature SNPs to reduce computation time
    """
    tmp_file = os.path.join(out_dir, "chr21_qsignature.vcf")
    if not utils.file_exists(tmp_file):
        cmd = ("grep chr21 {vcf_file} > {tmp_file}").format(**locals())
        out = subprocess.check_output(cmd, shell=True)
    return tmp_file
