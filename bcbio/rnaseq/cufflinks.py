"""Assess transcript abundance in RNA-seq experiments using Cufflinks.

http://cufflinks.cbcb.umd.edu/manual.html
"""
import os
import tempfile

from bcbio.utils import get_in, file_exists, safe_makedir
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
    options.extend(["--max-bundle-frags", 2000000])
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


def _get_output_dir(align_file, data, sample_dir=True):
    config = data["config"]
    name = data["rgnames"]["sample"] if sample_dir else ""
    return os.path.join(get_in(data, ("dirs", "work")), "cufflinks", name)

def assemble(bam_file, ref_file, gtf_file, num_cores, out_dir):
    safe_makedir(out_dir)
    with file_transaction(out_dir) as tmp_out_dir:
        cmd = ("cufflinks --output-dir {tmp_out_dir} --num-threads {num_cores} "
               "-g {gtf_file} --frag-bias-correct {ref_file} "
               "--multi-read-correct --upper-quartile-norm {bam_file}")
        cmd = cmd.format(**locals())
        do.run(cmd, "Assembling transcripts with Cufflinks using %s." % bam_file)
    return out_dir

def merge(assembled_gtfs, ref_file, gtf_file, num_cores):
#    assembled = " ".join(assembled_gtfs)
    assembled_file = tempfile.NamedTemporaryFile(delete=False).name
    with open(assembled_file, "w") as temp_handle:
        for assembled in assembled_gtfs:
            temp_handle.write(assembled + "\n")
    out_dir = os.path.join("assembly", "cuffmerge")
    out_file = os.path.join(out_dir, "merged.gtf")
    if file_exists(out_file):
        return out_file
    with file_transaction(out_dir) as tmp_out_dir:
        cmd = ("cuffmerge -o {tmp_out_dir} --ref-gtf {gtf_file} "
               "--num-threads {num_cores} --ref-sequence {ref_file} "
               "{assembled_file}")
        cmd = cmd.format(**locals())
        do.run(cmd, "Merging transcript assemblies with reference.")
    return out_file
