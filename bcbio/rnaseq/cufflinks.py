"""Assess transcript abundance in RNA-seq experiments using Cufflinks.

http://cufflinks.cbcb.umd.edu/manual.html
"""
import os
import subprocess

from bcbio.pipeline.variation import configured_ref_file

def assemble_transcripts(align_file, ref_file, config):
    """Create transcript assemblies using Cufflinks.
    """
    work_dir, fname = os.path.split(align_file)
    out_dir = os.path.join(work_dir,
                           "{base}-cufflinks".format(base=os.path.splitext(fname)[0]))
    cl = [config["program"].get("cufflinks", "cufflinks"),
          align_file,
          "-o", out_dir,
          "-b", ref_file,
          "-u"]
    tx_file = configured_ref_file("transcripts", config, ref_file)
    tx_mask_file = configured_ref_file("transcripts_mask", config, ref_file)
    if tx_file:
        cl += ["-g", tx_file]
    if tx_mask_file:
        cl += ["-M", tx_mask_file]
    out_tx_file = os.path.join(out_dir, "transcripts.gtf")
    if not os.path.exists(out_tx_file):
        subprocess.check_call(cl)
    assert os.path.exists(out_tx_file)
    return out_tx_file
