"""Filtering of genomic variants.
"""
from distutils.version import LooseVersion
import os

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import config_utils
from bcbio.provenance import do, programs

# ## General functionality

def jexl_hard(broad_runner, snp_file, ref_file, filter_type,
                                expressions):
    """Perform hard filtering with GATK using JEXL expressions.

    Variant quality score recalibration will not work on some regions; it
    requires enough positions to train the model. This provides a general wrapper
    around GATK to do cutoff based filtering.
    """
    base, ext = os.path.splitext(snp_file)
    out_file = "{base}-filter{ftype}{ext}".format(base=base, ext=ext,
                                                  ftype=filter_type)
    if not utils.file_exists(out_file):
        logger.debug("Hard filtering %s with %s" % (snp_file, expressions))
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "VariantFiltration",
                      "-R", ref_file,
                      "-l", "ERROR",
                      "--out", tx_out_file,
                      "--variant", snp_file]
            for exp in expressions:
                params.extend(["--filterName", "GATKStandard{e}".format(e=exp.split()[0]),
                               "--filterExpression", exp])
            broad_runner.run_gatk(params)
    return out_file

# ## Caller specific

def freebayes(in_file, ref_file, vrn_files, config):
    """FreeBayes filters: trying custom filter approach before falling back on hard filtering.
    """
    out_file = _freebayes_custom(in_file, ref_file, config)
    if out_file is None:
        out_file = _freebayes_hard(in_file, ref_file, config)
    return out_file

def _freebayes_custom(in_file, ref_file, config):
    """Custom FreeBayes filtering using bcbio.variation, tuned to human NA12878 results.
    """
    bv_ver = programs.get_version("bcbio.variation", config=config)
    if LooseVersion(bv_ver) < LooseVersion("0.1.1"):
        return None
    out_file = "%s-filter%s" % os.path.splitext(in_file)
    if not utils.file_exists(out_file):
        tmp_dir = utils.safe_makedir(os.path.join(os.path.dirname(in_file), "tmp"))
        bv_jar = config_utils.get_jar("bcbio.variation",
                                      config_utils.get_program("bcbio_variation", config, "dir"))
        resources = config_utils.get_resources("bcbio_variation", config)
        jvm_opts = resources.get("jvm_opts", ["-Xms750m", "-Xmx2g"])
        java_args = ["-Djava.io.tmpdir=%s" % tmp_dir]
        cmd = ["java"] + jvm_opts + java_args + ["-jar", bv_jar, "variant-filter", "freebayes",
                                                 in_file, ref_file]
        do.run(cmd, "Custom FreeBayes filtering using bcbio.variation")
    return out_file

def _freebayes_hard(in_file, ref_file, config):
    """Perform basic sanity filtering of FreeBayes results, removing low confidence calls.
    """
    broad_runner = broad.runner_from_config(config)
    filters = ["QUAL < 20.0", "DP < 5"]
    return jexl_hard(broad_runner, in_file, ref_file, "", filters)

