"""Call HLA alleles with assembly methods implemented in bwakit.

https://github.com/lh3/bwa/blob/master/README-alt.md#hla-typing
https://github.com/lh3/bwa/tree/master/bwakit
"""
import glob
import os

from bcbio import utils
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do

def run(data):
    bwakit_dir = os.path.dirname(os.path.realpath(utils.which("run-bwamem")))
    align_file = dd.get_align_bam(data)
    hla_base = os.path.join(utils.safe_makedir(os.path.join(os.path.dirname(align_file), "hla")),
                            os.path.basename(align_file) + ".hla")
    if len(glob.glob(hla_base + ".*")) > 0:
        out_file = hla_base + ".top"
        if not utils.file_exists(out_file):
            cmd = "{bwakit_dir}/run-HLA {hla_base}"
            do.run(cmd.format(**locals()), "HLA typing with bwakit")
        data["hla"] = {"top": out_file}
    return data
