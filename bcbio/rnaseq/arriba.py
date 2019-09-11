import os

from bcbio.heterogeneity import chromhacks
from bcbio.pipeline import config_utils, shared
from bcbio.pipeline import datadict as dd
from bcbio.log import logger
from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do

SUPPORTED_BUILDS = ("hg38", "GRCh37", "hg19")

def run_arriba(data):
    build = dd.get_genome_build(data)
    if build not in SUPPORTED_BUILDS:
        logger.info(f"{build} not supported for arriba, skipping.")
        return data

    arriba_dir = os.path.join(dd.get_work_dir(data), "arriba", dd.get_sample_name(data))
    utils.safe_makedir(arriba_dir)
    bam_file = dd.get_work_bam(data)
    ref_file = dd.get_ref_file(data)
    gtf = dd.get_gtf_file(data)
    arriba = config_utils.get_program("arriba", data)
    fusion_file = os.path.join(arriba_dir, "fusions.tsv")
    discarded_fusion_file = os.path.join(arriba_dir, "fusions.discarded.tsv")
    blacklist_file = get_arriba_blacklist_file(data)
    contigs = get_contigs(data)
    contig_list = ",".join(contigs)
    if utils.file_exists(fusion_file):
        data["arriba"] = {"fusions": fusion_file, "discarded": discarded_fusion_file}
        return(data)

    with file_transaction(fusion_file) as tx_fusion_file, \
         file_transaction(discarded_fusion_file) as tx_discarded_fusion_file:
        cmd = (f"{arriba} -x {bam_file} -g {gtf} -a {ref_file} -o {tx_fusion_file} "
               f"-O {tx_discarded_fusion_file} "
               f"-i {contig_list} ")
        if blacklist_file:
            logger.info(f"arriba blacklist file found, running blacklisting with {blacklist_file}.")
            cmd += (f"-b {blacklist_file} ")
        else:
            logger.info("arriba blacklist file not found, disabling blacklist filtering.")
            cmd += (f"-f blacklist ")
        message = f"Running arriba on {dd.get_sample_name(data)}."
        do.run(cmd, message)

    data["arriba"] = {"fusions": fusion_file, "discarded": discarded_fusion_file}
    return(data)

def get_arriba_blacklist_file(data):
    arriba_dir = os.path.join(os.path.dirname(dd.get_gtf_file(data)),
                              "fusion-blacklist")
    blacklist = os.path.join(arriba_dir, "arriba-blacklist.tsv.gz")
    if utils.file_exists(blacklist):
        return blacklist
    else:
        return None

def get_contigs(data):
    contigs = [x.name for x in shared.get_noalt_contigs(data)]
    keep = [x for x in contigs if chromhacks.is_autosomal(x) or chromhacks.is_sex(x)]
    return keep
