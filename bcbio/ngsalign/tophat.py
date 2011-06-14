"""Next-gen alignments with TopHat a spliced read mapper for RNA-seq experiments.

http://tophat.cbcb.umd.edu
"""
import os
import subprocess

galaxy_location_file = "bowtie_indices.loc"

def align(fastq_file, pair_file, ref_file, out_base, align_dir, config):
    out_dir = os.path.join(align_dir, "%s.sam" % out_base)
    if not os.path.exists(out_dir):
        cl = [config["program"]["tophat"]]
        cl += ["--solexa1.3-quals",
              "-p 8",
              "-r 45",
              ref_file]
        if pair_file:
            cl += [fastq_file, pair_file]
        else:
            cl += [fastq_file]
        cl += ["-o ", out_dir]
        child = subprocess.check_call(cl)
    return out_dir
