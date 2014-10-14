"""Assess transcript abundance in RNA-seq experiments using Cufflinks.

http://cufflinks.cbcb.umd.edu/manual.html
"""
import os

from bcbio.utils import get_in, file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.provenance import do


def run(align_file, ref_file, data):
    config = data["config"]
    cmd = _get_general_options(align_file, config)
    cmd.extend(_get_no_assembly_options(ref_file, data))
    out_dir = _get_output_dir(align_file, data)
    out_file = os.path.join(out_dir, "genes.fpkm_tracking")
    if file_exists(out_file):
        return out_dir
    with file_transaction(out_dir) as tmp_out_dir:
        cmd.extend(["--output-dir", tmp_out_dir])
        cmd.extend([align_file])
        cmd = map(str, cmd)
        do.run(cmd, "Cufflinks on %s." % (align_file))
    return out_dir

def _get_general_options(align_file, config):
    options = []
    cufflinks = config_utils.get_program("cufflinks", config)
    options.extend([cufflinks])
    options.extend(["--num-threads", config["algorithm"].get("num_cores", 1)])
    options.extend(["--quiet"])
    options.extend(["--no-update-check"])
    return options

def _get_no_assembly_options(ref_file, data):
    options = []
    options.extend(["--frag-bias-correct", ref_file])
    options.extend(["--multi-read-correct"])
    options.extend(["--upper-quartile-norm"])
    gtf_file = data["genome_resources"]["rnaseq"].get("transcripts", "")
    if gtf_file:
        options.extend(["--GTF", gtf_file])
    mask_file = data["genome_resources"]["rnaseq"].get("transcripts_mask", "")
    if mask_file:
        options.extend(["--mask-file", mask_file])

    return options


def _get_output_dir(align_file, data):
    config = data["config"]
    name = data["rgnames"]["sample"]
    return os.path.join(get_in(data, ("dirs", "work")), "cufflinks", name)
