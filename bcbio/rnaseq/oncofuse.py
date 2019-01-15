"""annotate fusion outputs from STAR and Tophat

Supported:
  oncofuse: http://www.unav.es/genetica/oncofuse.html
  github: https://github.com/mikessh/oncofuse
"""
from __future__ import print_function
import os
import pysam
import toolz as tz

from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.provenance import do
import bcbio.pipeline.datadict as dd
from bcbio.log import logger

# ## oncofuse fusion trancript detection

def run(data):
    if not aligner_supports_fusion(data):
        aligner = dd.get_aligner(data)
        logger.warning("Oncofuse is not supported for the %s aligner, "
                       "skipping. " % aligner)
        return None
    config = data["config"]
    genome_build = data.get("genome_build", "")
    input_type, input_dir, input_file = _get_input_para(data)
    if genome_build == "GRCh37":  # assume genome_build is hg19 otherwise
        if config["algorithm"].get("aligner") in ["star"]:
            if file_exists(input_file):
                input_file = _fix_star_junction_output(input_file)
        if config["algorithm"].get("aligner") in ["tophat", "tophat2"]:
            if file_exists(input_file):
                input_file = _fix_tophat_junction_output(input_file)
    elif "hg19" not in genome_build:
        return None
    #handle cases when fusion file doesn't exist
    if not file_exists(input_file):
        return None
    out_file = os.path.join(input_dir, "oncofuse_out.txt")
    if file_exists(out_file):
        return out_file
    oncofuse = config_utils.get_program("oncofuse", config)

    tissue_type = _oncofuse_tissue_arg_from_config(data)
    resources = config_utils.get_resources("oncofuse", config)
    if not file_exists(out_file):
        cl = [oncofuse]
        cl += resources.get("jvm_opts", ["-Xms750m", "-Xmx5g"])
        with file_transaction(data, out_file) as tx_out_file:
            cl += [input_file, input_type, tissue_type, tx_out_file]
            cmd = " ".join(cl)
            try:
                do.run(cmd, "oncofuse fusion detection", data)
            except:
                do.run("touch %s && echo '# failed' >> %s" % (tx_out_file, tx_out_file), "oncofuse failed", data)
                #return out_file
    return out_file

def is_non_zero_file(fpath):
    return True if os.path.isfile(fpath) and os.path.getsize(fpath) > 0 else False

def aligner_supports_fusion(data):
    SUPPORTED_ALIGNERS = ["tophat2", "tophat", "star"]
    aligner = dd.get_aligner(data)
    return aligner and aligner.lower() in SUPPORTED_ALIGNERS

def _get_input_para(data):
    TOPHAT_FUSION_OUTFILE = "fusions.out"
    STAR_FUSION_OUTFILE = "Chimeric.out.junction"
    config = data["config"]
    is_disambiguate = len(config["algorithm"].get("disambiguate", [])) > 0
    aligner = config["algorithm"].get("aligner")
    if aligner == "tophat2":
        aligner = "tophat"
    names = data["rgnames"]
    # set some default hard filters:
    N = 2 # min. spanning reads
    M = 4 # min. supporting reads (spanning + encompassing)
    align_dir_parts = os.path.join(data["dirs"]["work"], "align", names["sample"])
    align_dir_parts = os.path.join(align_dir_parts, data["genome_build"]) if is_disambiguate else align_dir_parts

    if aligner in ["tophat", "tophat2"]:
        align_dir_parts = os.path.join(data["dirs"]["work"], "align", names["sample"], names["lane"]+"_%s" % aligner)
        return "tophat-%d-%d" % (N,M), align_dir_parts, os.path.join(align_dir_parts, TOPHAT_FUSION_OUTFILE)
    if aligner in ["star"]:
        star_junction_file = os.path.join(align_dir_parts, names["lane"]+STAR_FUSION_OUTFILE)
        if is_disambiguate:
            contamination_bam = data["disambiguate"][ config["algorithm"]["disambiguate"][0] ]
            disambig_out_file = star_junction_file + "_disambiguated"
            if file_exists(disambig_out_file):
                star_junction_file = disambig_out_file
            elif file_exists(star_junction_file) and file_exists(contamination_bam):
                star_junction_file = _disambiguate_star_fusion_junctions(star_junction_file, contamination_bam,
                                                                         disambig_out_file, data)
        return "rnastar-%d-%d" % (N,M), align_dir_parts, star_junction_file
    return None

def _fix_tophat_junction_output(chimeric_out_junction_file):
    #for fusion.out
    out_file = chimeric_out_junction_file + ".hg19"
    with open(out_file, "w") as out_handle:
        with open(chimeric_out_junction_file, "r") as in_handle:
            for line in in_handle:
                parts = line.split("\t")
                left, right = parts[0].split("-")
                leftchr = _h37tohg19(left)
                rightchr = _h37tohg19(right)
                if not leftchr or not rightchr:
                    continue
                parts[0] = "%s-%s" % (_h37tohg19(left), _h37tohg19(right))
                out_handle.write("\t".join(parts))
    return out_file

def _fix_star_junction_output(chimeric_out_junction_file):
    #for Chimeric.out.junction
    out_file = chimeric_out_junction_file + ".hg19"
    with open(out_file, "w") as out_handle:
        with open(chimeric_out_junction_file, "r") as in_handle:
            for line in in_handle:
                parts = line.split("\t")
                parts[0] = _h37tohg19(parts[0])
                parts[3] = _h37tohg19(parts[3])
                if not parts[0] or not parts[3]:
                    continue
                out_handle.write("\t".join(parts))
    return out_file

def _h37tohg19(chromosome):
    MAX_CHROMOSOMES = 23
    if chromosome in [str(x) for x in range(1, MAX_CHROMOSOMES)] + ["X", "Y"]:
        new_chrom = "chr%s" % chromosome
    elif chromosome == "MT":
        new_chrom = "chrM"
    # not a supported chromosome
    else:
        return None
    return new_chrom


def _oncofuse_tissue_arg_from_config(data):

    """Retrieve oncofuse arguments supplied through input configuration.
    tissue_type is the library argument, which tells Oncofuse to use its
    own pre-built gene expression libraries. There are four pre-built
    libraries, corresponding to the four supported tissue types:
    EPI (epithelial origin),
    HEM (hematological origin),
    MES (mesenchymal origin) and
    AVG (average expression, if tissue source is unknown).
    """
    SUPPORTED_TISSUE_TYPE = ["EPI", "HEM", "MES", "AVG"]
    tissue_type = tz.get_in(("metadata", "tissue"), data, None)
    if not tissue_type:
        logger.info("Oncofuse: tissue type not set, using average expression (AVG).")
        tissue_type = "AVG"
    elif tissue_type not in SUPPORTED_TISSUE_TYPE:
        logger.info("Oncofuse: %s not a supported tissue type, using average "
                    "expression (AVG)." % tissue_type)
        tissue_type = "AVG"
    else:
        logger.info("Oncofuse: using %s as tissue type." % tissue_type)
    return tissue_type

def _disambiguate_star_fusion_junctions(star_junction_file, contamination_bam, disambig_out_file, data):
    """ Disambiguate detected fusions based on alignments to another species.
    """
    out_file = disambig_out_file
    fusiondict = {}
    with open(star_junction_file, "r") as in_handle:
        for my_line in in_handle:
            my_line_split = my_line.strip().split("\t")
            if len(my_line_split) < 10:
                continue
            fusiondict[my_line_split[9]] = my_line.strip("\n")
    with pysam.Samfile(contamination_bam, "rb") as samfile:
        for my_read in samfile:
            if my_read.is_unmapped or my_read.is_secondary:
                continue
            if my_read.qname in fusiondict:
                fusiondict.pop(my_read.qname)
    with file_transaction(data, out_file) as tx_out_file:
        with open(tx_out_file, 'w') as myhandle:
            for my_key in fusiondict:
                print(fusiondict[my_key], file=myhandle)

    return out_file
