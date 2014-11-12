"""identify somatic variants in cancer samples

https://qcmg.org/bioinformatics/tiki-index.php?page=qSNP#EXAMPLES
"""

from __future__ import print_function
import os
import shutil
from re import sub

try:
    import vcf
except ImportError:
    vcf = None

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import subset_variant_regions, remove_lcr_regions
from bcbio.provenance import do
from bcbio.variation import annotation
from bcbio.variation.vcfutils import get_paired_bams, is_paired_analysis, bgzip_and_index

def is_installed(config):
    """Check for qsnp installation on machine.
    """
    try:
        config_utils.get_program("qsnp", config)
        return True
    except config_utils.CmdNotFound:
        return False


def run_qsnp(align_bams, items, ref_file, assoc_files, region=None,
                  out_file=None):
    """Run Qsnp indel calling, either paired tumor/normal or germline calling.
    """
    if is_paired_analysis(align_bams, items):
        call_file = _run_qsnp_paired(align_bams, items, ref_file,
                                          assoc_files, region, out_file)
        return call_file
    return None


def _run_qsnp_paired(align_bams, items, ref_file, assoc_files,
                          region=None, out_file=None):
    """Detect indels with Qsnp.

    This is used for paired tumor / normal samples.
    """
    config = items[0]["config"]
    if out_file is None:
        out_file = "%s-paired-variants.vcf" % os.path.splitext(align_bams[0])[0]
    if not utils.file_exists(out_file):
        out_file = out_file.replace(".gz", "")
        with file_transaction(config, out_file) as tx_out_file:
            with tx_tmpdir() as tmpdir:
                paired = get_paired_bams(align_bams, items)
                if not paired.normal_bam:
                    return None
                qsnp = config_utils.get_program("qsnp", config)
                qsnp_log = os.path.join(tmpdir, "qsnp.log")
                qsnp_init = os.path.join(tmpdir, "qsnp.ini")
                _create_input(paired, tx_out_file, ref_file, assoc_files['dbsnp'], qsnp_init)

                cl = ("{qsnp} -i {qsnp_init} -log {qsnp_log}")
                do.run(cl.format(**locals()), "Genotyping paired variants with Qsnp", {})
        out_file = _fix_vcf(out_file)
        out_file = bgzip_and_index(out_file, config)
    return out_file


def _create_input(paired, out_file, ref_file, snp_file, qsnp_file):
    ini_file["[inputFiles]"]["dbSNP"] = snp_file
    ini_file["[inputFiles]"]["ref"] = ref_file
    ini_file["[inputFiles]"]["normalBam"] = paired.normal_bam
    ini_file["[inputFiles]"]["tumourBam"] = paired.tumor_bam
    ini_file["[ids]"]["normalSample"] = paired.normal_name
    ini_file["[ids]"]["tumourSample"] = paired.tumor_name
    ini_file["[ids]"]["donor"] = paired.tumor_name
    ini_file["[outputFiles]"]["vcf"] = out_file
    with open(qsnp_file, "w") as out_handle:
        for k, v in ini_file.iteritems():
            out_handle.write("%s\n" % k)
            for opt, value in v.iteritems():
                if value != "":
                    out_handle.write("%s = %s\n" % (opt, value))


def _fix_vcf(out_file):
    """Fix FILTER and FORMAT fields"""
    in_file = out_file.replace(".vcf", "-ori.vcf")
    SBIAS_line = ('##FILTER=<ID=SBIAS,Description="Due to bias">\n'
                  '##FILTER=<ID=5BP,Description="Due to 5BP">\n')
    SOMATIC_line = '##INFO=<ID=SOMATIC,Number=1,Type=String,Description="somatic variation">\n'
    if not utils.file_exists(in_file):
        shutil.move(out_file, in_file)
    with file_transaction(out_file) as tx_out_file:
        with open(in_file) as in_handle, open(tx_out_file, "w") as out_handle:
            for line in in_handle:
                if "5BP" in line:
                    line = sub("5BP[0-9]+", "5BP", line)
                out_handle.write("%s" % line)
                if SBIAS_line and line.startswith("##FILTER"):
                    out_handle.write("%s" % SBIAS_line)
                    SBIAS_line = ""
                if line.startswith("##INFO") and SOMATIC_line:
                    out_handle.write("%s" % SOMATIC_line)
                    SOMATIC_line = ""
    return out_file


ini_file= {"[inputFiles]":{
          "dbSNP":"",
          "germlineDB":"",
          "chrConv":"",
          "ref":"",
          "normalBam":"",
          "tumourBam":""},

          "[parameters]":{
          "annotateMode":"vcf",
          "runMode":"standard",
          "minimumBaseQuality":"10",
          "includeIndels": "true"},

          "[ids]":{
          "donor":"",
          "normalSample":"",
          "tumourSample":"",
          "somaticAnalysis":"",
          "germlineAnalysis":""},

          "[outputFiles]":{
          "vcf":"",
          "dccSomatic":"",
          "dccGermline":""},

          "[rules]":{
          "normal1":"0,20,3",
          "normal2":"21,50,4",
          "normal3":"51,,10",
          "tumour1":"0,20,3",
          "tumour2":"21,50,4",
          "tumour3":"51,,5"}}



