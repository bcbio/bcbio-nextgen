import os
from bcbio import utils
from bcbio.log import logger
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.pipeline import config_utils
import bcbio.pipeline.datadict as dd

def vcfanno(vcf, out_file, conf_fns, data, basepath=None, lua_fns=None):
    """
    annotate a VCF file using vcfanno (https://github.com/brentp/vcfanno)
    """
    if utils.file_exists(out_file):
        return out_file
    if lua_fns is None:
        lua_fns = []
    vcfanno = config_utils.get_program("vcfanno", data)
    with file_transaction(out_file) as tx_out_file:
        conffn = _combine_files(conf_fns, tx_out_file)
        luafn = _combine_files(lua_fns, tx_out_file)
        luaflag = "-lua {0}".format(luafn) if luafn and utils.file_exists(luafn) else ""
        basepathflag = "-base-path {0}".format(basepath) if basepath else ""
        cores = dd.get_num_cores(data)
        cmd = "{vcfanno} -p {cores} {luaflag} {basepathflag} {conffn} {vcf} | bgzip -c > {tx_out_file}"
        message = "Annotating {vcf} with vcfanno, using {conffn}".format(**locals())
        do.run(cmd.format(**locals()), message)
    return out_file

def _combine_files(orig_files, base_out_file):
    orig_files = [x for x in orig_files if x and utils.file_exists(x)]
    if not orig_files:
        return None
    elif len(orig_files) == 1:
        return orig_files[0]
    else:
        out_file = "%s-combine%s" % (utils.splitext_plus(base_out_file)[0],
                                     utils.splitext_plus(orig_files[0])[-1])
        with open(out_file, "w") as out_handle:
            for orig_file in orig_files:
                with open(orig_file) as in_handle:
                    for line in in_handle:
                        out_handle.write(line)
                out_handle.write("\n\n")
        return out_file

def run_vcfanno(vcf, conf_files, data, data_basepath=None):
    """
    annotated a VCF file using vcfanno, looks up the proper config/lua scripts
    under the `vcfanno` key under the algorithm section of the datadict,
    skipping if the files cannot be found
    """
    if not isinstance(conf_files, (list, tuple)):
        conf_files = [conf_files]
    build = dd.get_genome_build(data)
    basepath = os.path.abspath(os.path.join(os.path.dirname(dd.get_ref_file(data)),
                                            os.pardir))
    annodir = os.path.abspath(os.path.join(basepath, "config", "vcfanno"))
    conf_fns = []
    lua_fns = []
    anno_type = None
    for conf_file in conf_files:
        if utils.file_exists(conf_file) and os.path.isfile(conf_file):
            conffn = conf_file
            luafn = "%s.lua" % utils.splitext_plus(conffn)[0]
        else:
            anno_type = os.path.basename(conf_file)
            conffn = os.path.join(annodir, anno_type + ".conf")
            luafn = os.path.join(annodir, anno_type + ".lua")
        if not utils.file_exists(conffn):
            CONF_NOT_FOUND = (
                "The vcfanno configuration {conffn} was not found for {build}, skipping.")
            logger.warn(CONF_NOT_FOUND.format(**locals()))
        else:
            conf_fns.append(conffn)
            lua_fns.append(luafn)
    if not conf_fns:
        return vcf
    if not anno_type:
        anno_type = "gemini"
    out_file = utils.splitext_plus(vcf)[0] + "-annotated-" + anno_type + ".vcf.gz"
    if utils.file_exists(out_file):
        return out_file

    out_file = vcfanno(vcf, out_file, conf_fns, data, data_basepath or basepath, lua_fns)
    return out_file
