"""Next-gen sequencing alignment with Novoalign: http://www.novocraft.com

For BAM input handling this requires:
  novoalign (with license for multicore)
  novosort (also with license for multicore)
  samtools
"""
import os
import subprocess

from bcbio.pipeline import config_utils
from bcbio.log import logger
from bcbio.utils import (memoize_outfile, file_exists, transform_to, curdir_tmpdir)
from bcbio.distributed.transaction import file_transaction

# ## BAM realignment

def get_rg_info(names):
    return r"@RG\tID:{rg}\tPL:{pl}\tPU:{pu}\tSM:{sample}".format(**names)

def align_bam(in_bam, ref_file, names, align_dir, config):
    """Perform realignment of input BAM file, handling sorting of input/output with novosort.

    Uses unix pipes for avoid IO writing between steps:
      - novosort of input BAM to coordinates
      - alignment with novoalign
      - conversion to BAM with samtools
      - coordinate sorting with novosort
    """
    out_file = os.path.join(align_dir, "{0}-sort.bam".format(names["lane"]))
    novosort = config_utils.get_program("novosort", config)
    novoalign = config_utils.get_program("novoalign", config)
    samtools = config_utils.get_program("samtools", config)
    resources = config_utils.get_resources("novoalign", config)
    num_cores = config["algorithm"].get("num_cores", 1)
    max_mem = resources.get("memory", "4G")

    if not file_exists(out_file):
        with curdir_tmpdir(base_dir=align_dir) as work_dir:
            with file_transaction(out_file) as tx_out_file:
                rg_info = get_rg_info(names)
                cmd = ("{novosort} -c {num_cores} -m {max_mem} --compression 0 "
                       " -n -t {work_dir} {in_bam} "
                       "| {novoalign} -o SAM '{rg_info}' -d {ref_file} -f /dev/stdin "
                       "  -F BAMPE -c {num_cores} "
                       "| {samtools} view -b -S -u - "
                       "| {novosort} -c {num_cores} -m {max_mem} -t {work_dir} "
                       "  -o {tx_out_file} /dev/stdin")
                subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

# ## Fastq to BAM alignment

def _novoalign_args_from_config(config):
    """Select novoalign options based on configuration parameters.
    """
    qual_format = config["algorithm"].get("quality_format", "").lower()
    qual_flags = ["-F", "ILMFQ" if qual_format == "illumina" else "STDFQ"]
    multi_mappers = config["algorithm"].get("multiple_mappers", True)
    if multi_mappers is True:
        multi_flag = "Random"
    elif isinstance(multi_mappers, basestring):
        multi_flag = multi_mappers
    else:
        multi_flag = "None"
    multi_flags = ["-r"] + multi_flag.split()
    extra_args = config["algorithm"].get("extra_align_args", [])
    return qual_flags + multi_flags + extra_args

# Tweaks to add
# -k -t 200 -K quality calibration metrics
# paired end sizes

def align(fastq_file, pair_file, ref_file, out_base, align_dir, config,
          extra_args=None, rg_name=None):
    """Align with novoalign.
    """
    out_file = os.path.join(align_dir, "{0}.sam".format(out_base))
    if not file_exists(out_file):
        cl = [config_utils.get_program("novoalign", config)]
        cl += _novoalign_args_from_config(config)
        cl += extra_args if extra_args is not None else []
        cl += ["-o", "SAM"]
        if rg_name:
            cl.append(r"@RG\tID:{0}".format(rg_name))
        cl += ["-d", ref_file, "-f", fastq_file]
        if pair_file:
            cl.append(pair_file)
        with file_transaction(out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                logger.info(" ".join([str(x) for x in cl]))
                subprocess.check_call([str(x) for x in cl], stdout=out_handle)
    return out_file

# ## Indexing

@memoize_outfile(ext=".ndx")
def refindex(ref_file, kmer_size=None, step_size=None, out_file=None):
    cl = ["novoindex"]
    if kmer_size:
        cl += ["-k", str(kmer_size)]
    if step_size:
        cl += ["-s", str(step_size)]
    cl += [out_file, ref_file]
    subprocess.check_call(cl)

def remap_index_fn(ref_file):
    """Map sequence references to equivalent novoalign indexes.
    """
    checks = [os.path.splitext(ref_file)[0].replace("/seq/", "/novoalign/"),
              os.path.splitext(ref_file)[0] + ".ndx",
              ref_file + ".bs.ndx",
              ref_file + ".ndx"]
    for check in checks:
        if os.path.exists(check):
            return check
    return checks[0]
