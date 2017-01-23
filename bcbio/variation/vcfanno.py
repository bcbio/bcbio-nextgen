import os
from bcbio import utils
from bcbio.log import logger
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.pipeline import config_utils
import bcbio.pipeline.datadict as dd

def vcfanno(vcf, out_file, conffn, data, basepath=None, luafn=None):
    """
    annotate a VCF file using vcfanno (https://github.com/brentp/vcfanno)
    """
    if utils.file_exists(out_file):
        return out_file
    vcfanno = config_utils.get_program("vcfanno", data)
    luaflag = "-lua {0}".format(luafn) if luafn and utils.file_exists(luafn) else ""
    basepathflag = "-base-path {0}".format(basepath) if basepath else ""
    cores = dd.get_num_cores(data)
    cmd = "{vcfanno} -p {cores} {luaflag} {basepathflag} {conffn} {vcf} | bgzip -c > {tx_out_file}"
    message = "Annotating {vcf} with vcfanno, using {conffn}".format(**locals())
    with file_transaction(out_file) as tx_out_file:
        do.run(cmd.format(**locals()), message)
    return out_file

def run_vcfanno(vcf, anno_type, data, data_basepath=None):
    """
    annotated a VCF file using vcfanno, looks up the proper config/lua scripts
    under the `vcfanno` key under the algorithm section of the datadict,
    skipping if the files cannot be found
    """
    build = dd.get_genome_build(data)
    basepath = os.path.abspath(os.path.join(os.path.dirname(dd.get_ref_file(data)),
                                            os.pardir))
    annodir = os.path.abspath(os.path.join(basepath, "config", "vcfanno"))
    conffn = os.path.join(annodir, anno_type + ".conf")
    luafn = os.path.join(annodir, anno_type + ".lua")
    CONF_NOT_FOUND = (
        "The vcfanno configuration {conffn} was not found for {build}, skipping.")
    if not utils.file_exists(conffn):
        logger.warn(CONF_NOT_FOUND.format(**locals()))
        return vcf

    out_file = utils.splitext_plus(vcf)[0] + "-annotated-" + anno_type + ".vcf.gz"
    if utils.file_exists(out_file):
        return out_file

    out_file = vcfanno(vcf, out_file, conffn, data, data_basepath or basepath, luafn)
    return out_file
