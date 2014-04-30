"""annonate fusion transcript using external programs.

Supported:
  oncofuse: http://www.unav.es/genetica/oncofuse.html
"""

import os
import csv
import glob

from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.provenance import do

# ## oncofuse fusion trancript detection
#haven't tested with STAR, instructions referenced from seqanswer, http://seqanswers.com/forums/archive/index.php/t-33095.html

def run(data):
    #cmd line: java -Xmx1G -jar Oncofuse.jar input_file input_type tissue_type output_file
    config = data["config"]
    genome_build = data.get("genome_build", "")
    input_type, input_dir, input_file = _get_input_para(data)
    if genome_build == 'GRCh37': #assume genome_build is hg19 otherwise
        if config["algorithm"].get("aligner") in ['star']:
            input_file = _fix_star_junction_output(input_file)
        if config["algorithm"].get("aligner") in ['tophat', 'tophat2']:
            input_file = _fix_tophat_junction_output(input_file)
    
    #handle cases when fusion file doesn't exist
    if not file_exists(input_file):
        return None
    
    out_file = os.path.join(input_dir, 'oncofuse_out.txt')
    
    if file_exists(out_file):
        return out_file
    
    oncofuse_jar = config_utils.get_jar("Oncofuse",
                                      config_utils.get_program("oncofuse",
                                                               config, "dir"))

    tissue_type = _oncofuse_tissue_arg_from_config(data)
    resources = config_utils.get_resources("oncofuse", config)
    if not file_exists(out_file):
        cl = ["java"]
        cl += resources.get("jvm_opts", ["-Xms750m", "-Xmx5g"])
        cl += ["-jar", oncofuse_jar, input_file, input_type, tissue_type, out_file]
        with open(out_file, "w") as out_handle:
            cmd = " ".join(cl)
            try:
                do.run(cmd, "oncofuse fusion detection", data)
            except:
                return out_file
    return out_file

def is_non_zero_file(fpath):  
    return True if os.path.isfile(fpath) and os.path.getsize(fpath) > 0 else False

def _get_input_para(data):

    TOPHAT_FUSION_OUTFILE = "fusions.out"
    STAR_FUSION_OUTFILE = 'Chimeric.out.junction'
    
    
    config = data["config"]
    aligner = config["algorithm"].get("aligner")
    if aligner == 'tophat2':
        aligner = 'tophat'
    names = data["rgnames"]
    align_dir_parts = os.path.join(data["dirs"]["work"], "align", names["lane"], names["sample"]+"_%s" % aligner)
    if aligner in ['tophat', 'tophat2']:
        align_dir_parts = os.path.join(data["dirs"]["work"], "align", names["lane"], names["sample"]+"_%s" % aligner)
        return 'tophat', align_dir_parts, os.path.join(align_dir_parts, TOPHAT_FUSION_OUTFILE)
    if aligner in ['star']:
        align_dir_parts = os.path.join(data["dirs"]["work"], "align", names["lane"])
        return 'rnastar', align_dir_parts, os.path.join(align_dir_parts,names["lane"]+STAR_FUSION_OUTFILE)
    return None

def _fix_tophat_junction_output(chimeric_out_junction_file):
    #for fusion.out
    out_file = chimeric_out_junction_file + '.hg19'
    with open(out_file, "w") as out_handle:
        with open(chimeric_out_junction_file, "r") as in_handle:
            for line in in_handle:
                parts = line.split("\t")
                left, right = parts[0].split("-")
                parts[0] = "%s-%s" % (_h37tohg19(left), _h37tohg19(right))
                out_handle.write("\t".join(parts))
    return out_file    
    
def _fix_star_junction_output(chimeric_out_junction_file):
    #for Chimeric.out.junction
    out_file = chimeric_out_junction_file + '.hg19'
    with open(out_file, "w") as out_handle:
        with open(chimeric_out_junction_file, "r") as in_handle:
            for line in in_handle:
                parts = line.split("\t")
                parts[0] = _h37tohg19(parts[0])
                parts[3] = _h37tohg19(parts[3])
                out_handle.write("\t".join(parts))
    return out_file

def _h37tohg19(chromosome):
    MAX_CHROMOSOMES = 23
    if chromosome in [str(x) for x in range(1, MAX_CHROMOSOMES)] + ["X", "Y"]:
        new_chrom = "chr%s" % chromosome
    elif chromosome == "MT":
        new_chrom = "chrM"
    else:
        raise NotImplementedError(chromosome)
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
    SUPPORTED_TIISUE_TYPE = ["EPI", "HEM", "MES", "AVG"]
    if data.get("metadata", {}).get("tissue") in SUPPORTED_TIISUE_TYPE:
        return data.get("metadata", {}).get("tissue")
    else:
        return 'AVG'