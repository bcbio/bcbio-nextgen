"""annonate fusion transcript using external programs.

Supported:
  oncofuse: http://www.unav.es/genetica/oncofuse.html
"""

import os
import csv
import glob
import subprocess

from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils


# ## oncofuse fusion trancript detection
#haven't tested with STAR, instructions referenced from seqanswer, http://seqanswers.com/forums/archive/index.php/t-33095.html

def run(data):

    config = data["config"]
    input_type, input_dir, input_file = _get_input_para(data)
    out_file = os.path.join(input_dir, 'oncofuse_out.txt')
    #cmd line: java -Xmx1G -jar Oncofuse.jar input_file input_type tissue_type output_file
    oncofuse_jar = config_utils.get_jar("oncofuse",
                                      config_utils.get_program("oncofuse", config, "dir"))
    
    tissue_type = _oncofuse_tissue_arg_from_config(data)
    resources = config_utils.get_resources("oncofuse", config)
    if not file_exists(out_file):
        cl = ["java"]
        cl += resources.get("jvm_opts", ["-Xms750m", "-Xmx5g"])
        cl += ["-jar", oncofuse_jar, input_file, input_type, tissue_type, out_file]
        with file_transaction(out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                subprocess.check_call(cl, stdout=out_handle)
    return out_file


def _get_input_para(data):
       
    TOPHAT_FUSION_OUTFILE = "fusions.out"
    STAR_FUSION_OUTFILE = 'Chimeric.out.junction.txt'
    
    config = data["config"]
    aligner = config["algorithm"].get("aligner")
    align_dir_parts = [data["dirs"]["work"], "align", names["sample"]]
    
    if aligner in ['tophat', 'tophat2']:
        return 'tophat', align_dir_parts, os.path.join(align_dir_parts, TOPHAT_FUSION_OUTFILE)
    if aligner in ['star']:
        return 'rnastar', align_dir_parts, os.path.join(align_dir_parts, STAR_FUSION_OUTFILE)
    return None 



def _oncofuse_tissue_arg_from_config(data):
    
    """Retrieve oncofuse arguments supplied through input configuration.
    tissue_type is the library argument, which tells Oncofuse to use its own pre-built gene expression libraries.
    There are four pre-built libraries, corresponding to the four supported tissue types:
    EPI (epithelial origin),
    HEM (hematological origin),
    MES (mesenchymal origin) and
    AVG (average expression, if tissue source is unknown).
    """
    #potential check for tumor only analysis
    #if data.get("metadata", {}).get("tissue") in ["tumor", "normal"]:
    #    pass
    SUPPORTED_TIISUE_TYPE = ["EPI", "HEM", "MES", "AVG"]
    if data.get("metadata", {}).get("tissue") in SUPPORTED_TIISUE_TYPE:
        return [data.get("metadata", {}).get("tissue")]
    else:
        #may handle exception later
        return ['AVG']



