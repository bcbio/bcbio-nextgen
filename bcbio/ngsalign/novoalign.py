"""Next-gen sequencing alignment with Novoalign: http://www.novocraft.com

For BAM input handling this requires:
  novoalign (with license for multicore)
  samtools
"""
import os
import subprocess

from bcbio import bam, utils
from bcbio.ngsalign import alignprep, postalign
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.distributed.transaction import tx_tmpdir
from bcbio.utils import (memoize_outfile, file_exists)

# ## BAM realignment

def get_rg_info(names):
    return r"@RG\tID:{rg}\tPL:{pl}\tPU:{pu}\tSM:{sample}".format(**names)

def align_bam(in_bam, ref_file, names, align_dir, data):
    """Perform realignment of input BAM file; uses unix pipes for avoid IO.
    """
    config = data["config"]
    out_file = os.path.join(align_dir, "{0}-sort.bam".format(names["lane"]))
    novoalign = config_utils.get_program("novoalign", config)
    samtools = config_utils.get_program("samtools", config)
    resources = config_utils.get_resources("novoalign", config)
    num_cores = config["algorithm"].get("num_cores", 1)
    max_mem = resources.get("memory", "4G").upper()
    extra_novo_args = " ".join(_novoalign_args_from_config(config, False))

    if not file_exists(out_file):
        with tx_tmpdir(data, base_dir=align_dir) as work_dir:
            with postalign.tobam_cl(data, out_file, bam.is_paired(in_bam)) as (tobam_cl, tx_out_file):
                rg_info = get_rg_info(names)
                tx_out_prefix = os.path.splitext(tx_out_file)[0]
                prefix1 = "%s-in1" % tx_out_prefix
                cmd = ("{samtools} sort -n -o -l 1 -@ {num_cores} -m {max_mem} {in_bam} {prefix1} "
                       "| {novoalign} -o SAM '{rg_info}' -d {ref_file} -f /dev/stdin "
                       "  -F BAMPE -c {num_cores} {extra_novo_args} | ")
                cmd = (cmd + tobam_cl).format(**locals())
                do.run(cmd, "Novoalign: %s" % names["sample"], None,
                       [do.file_nonempty(tx_out_file), do.file_reasonable_size(tx_out_file, in_bam)])
    return out_file

# ## Fastq to BAM alignment

def align_pipe(fastq_file, pair_file, ref_file, names, align_dir, data):
    """Perform piped alignment of fastq input files, generating sorted output BAM.
    """
    pair_file = pair_file if pair_file else ""
    out_file = os.path.join(align_dir, "{0}-sort.bam".format(names["lane"]))
    if data.get("align_split"):
        final_file = out_file
        out_file, data = alignprep.setup_combine(final_file, data)
        fastq_file = alignprep.split_namedpipe_cl(fastq_file, data)
        if pair_file:
            pair_file = alignprep.split_namedpipe_cl(pair_file, data)
    else:
        final_file = None
    samtools = config_utils.get_program("samtools", data["config"])
    novoalign = config_utils.get_program("novoalign", data["config"])
    resources = config_utils.get_resources("novoalign", data["config"])
    num_cores = data["config"]["algorithm"].get("num_cores", 1)
    max_mem = resources.get("memory", "1G")
    extra_novo_args = " ".join(_novoalign_args_from_config(data["config"]))
    rg_info = get_rg_info(names)
    if not utils.file_exists(out_file) and (final_file is None or not utils.file_exists(final_file)):
        with tx_tmpdir(data) as work_dir:
            with postalign.tobam_cl(data, out_file, pair_file != "") as (tobam_cl, tx_out_file):
                tx_out_prefix = os.path.splitext(tx_out_file)[0]
                cmd = ("{novoalign} -o SAM '{rg_info}' -d {ref_file} -f {fastq_file} {pair_file} "
                       "  -c {num_cores} {extra_novo_args} | ")
                cmd = (cmd + tobam_cl).format(**locals())
                do.run(cmd, "Novoalign: %s" % names["sample"], None,
                       [do.file_nonempty(tx_out_file), do.file_reasonable_size(tx_out_file, fastq_file)])
    data["work_bam"] = out_file
    return data

def _novoalign_args_from_config(config, need_quality=True):
    """Select novoalign options based on configuration parameters.
    """
    if need_quality:
        qual_format = config["algorithm"].get("quality_format", "").lower()
        qual_flags = ["-F", "ILMFQ" if qual_format == "illumina" else "STDFQ"]
    else:
        qual_flags = []
    multi_mappers = config["algorithm"].get("multiple_mappers")
    if multi_mappers is True:
        multi_flag = "Random"
    elif isinstance(multi_mappers, basestring):
        multi_flag = multi_mappers
    else:
        multi_flag = "None"
    multi_flags = ["-r"] + multi_flag.split()
    resources = config_utils.get_resources("novoalign", config)
    # default arguments for improved variant calling based on
    # comparisons to reference materials: turn off soft clipping and recalibrate
    if resources.get("options") is None:
        extra_args = ["-o", "FullNW", "-k"]
    else:
        extra_args = [str(x) for x in resources.get("options", [])]
    return qual_flags + multi_flags + extra_args

# Tweaks to add
# -k -t 200 -K quality calibration metrics
# paired end sizes

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

# ## Galaxy integration

# Optional galaxy location file. Falls back on remap_index_fn if not found
galaxy_location_file = "novoalign_indices.loc"

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
