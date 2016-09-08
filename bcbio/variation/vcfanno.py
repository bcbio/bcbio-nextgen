import os
from bcbio import utils
from bcbio.log import logger
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.pipeline import config_utils
import bcbio.pipeline.datadict as dd

SUPPORTED_ANNOTATION_TYPES = ["rnaedit"]

def vcfanno(vcf, out_file, conffn, data, basepath=None, luafn=None):
    """
    annotate a VCF file using vcfanno (https://github.com/brentp/vcfanno)
    """
    if utils.file_exists(out_file):
        return out_file
    vcfanno = config_utils.get_program("vcfanno", data)
    luaflag = "-lua {0}".format(luafn) if luafn else ""
    basepathflag = "-base-path {0}".format(basepath) if basepath else ""
    cmd = "{vcfanno} {luaflag} {basepathflag} {conffn} {vcf} | bgzip -c > {tx_out_file}"
    message = "Annotating {vcf} with vcfanno, using {conffn}".format(**locals())
    with file_transaction(out_file) as tx_out_file:
        do.run(cmd.format(**locals()), message)
    return out_file

def run_vcfanno(vcf, anno_type, data):
    """
    annotated a VCF file using vcfanno, looks up the proper config/lua scripts
    under the `vcfanno` key under the algorithm section of the datadict,
    skipping if the files cannot be found
    """
    UNSUPPORTED_TYPE_MESSAGE = (
        "{anno_type} is not a supported vcf annotation type with vcfanno. "
        "Supported types are {SUPPORTED_ANNOTATION_TYPES}")
    if anno_type not in SUPPORTED_ANNOTATION_TYPES:
        logger.warn(UNSUPPORTED_TYPE_MESSAGE.format(**locals()))
        return vcf
    build = dd.get_genome_build(data)
    annodir = os.path.dirname(dd.get_ref_file(data))
    annodir = os.path.abspath(os.path.join(annodir, os.pardir, "vcfanno"))
    annostem = os.path.join(annodir, build + "-")
    conffn = annostem + anno_type + ".conf"
    luafn = annostem + anno_type + ".lua"
    CONF_NOT_FOUND = (
        "The vcfanno configuration {conffn} was not found for {build}, skipping.")
    if not utils.file_exists(conffn):
        logger.warn(CONF_NOT_FOUND.format(**locals()))
        return vcf

    base = os.path.splitext(vcf)[0]
    out_file = base + anno_type + "-annotated.vcf.gz"
    if utils.file_exists(out_file):
        return out_file
    basepath = os.path.abspath(os.path.join(os.path.dirname(dd.get_ref_file(data)),
                                            os.path.pardir))
    basepath = annodir

    out_file = vcfanno(vcf, out_file, conffn, data, basepath, luafn)
    return out_file
