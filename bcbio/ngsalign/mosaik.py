"""Next gen sequence alignment with Mosaik.

http://bioinformatics.bc.edu/marthlab/Mosaik
"""
import os
import subprocess

from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction

galaxy_location_file = "mosaik_index.loc"

def _mosaik_args_from_config(config):
    """Configurable high level options for mosaik.
    """
    multi_mappers = config["algorithm"].get("multiple_mappers", True)
    multi_flags = ["-m", "all"] if multi_mappers else ["-m", "unique"]
    error_flags = ["-mm", config["algorithm"]["max_errors"]]
    cores = config.get("resources", {}).get("mosaik", {}).get("cores", None)
    core_flags = ["-p", str(cores)] if cores else []
    return core_flags + multi_flags + error_flags

def _convert_fastq(fastq_file, pair_file, rg_name, out_file, config):
    """Convert fastq inputs into internal Mosaik representation.
    """
    out_file = "{0}-fq.mkb".format(os.path.splitext(out_file)[0])
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            cl = [config["program"].get("mosaik", "MosaikAligner").replace("Aligner", "Build")]
            cl += ["-q", fastq_file,
                   "-out", tx_out_file,
                   "-st", config["algorithm"].get("platform", "illumina").lower()]
            if pair_file:
                cl += ["-q2", pair_file]
            if rg_name:
                cl += ["-id", rg_name]
            env_set = "export MOSAIK_TMP={0}".format(os.path.dirname(tx_out_file))
            subprocess.check_call(env_set + " && " + " ".join(cl), shell=True)
    return out_file

def _get_mosaik_nn_args(out_file):
    """Retrieve default neural network files from GitHub to pass to Mosaik.
    """
    base_nn_url = "https://raw.github.com/wanpinglee/MOSAIK/master/src/networkFile/"
    out = []
    for arg, fname in [("-annse", "2.1.26.se.100.005.ann"),
                       ("-annpe", "2.1.26.pe.100.0065.ann")]:
        arg_fname = os.path.join(os.path.dirname(out_file), fname)
        if not file_exists(arg_fname):
            subprocess.check_call(["wget", "-O", arg_fname, base_nn_url + fname])
        out += [arg, arg_fname]
    return out

def align(fastq_file, pair_file, ref_file, out_base, align_dir, config,
          extra_args=None, rg_name=None):
    """Alignment with MosaikAligner.
    """
    out_file = os.path.join(align_dir, "%s-align.bam" % out_base)
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            built_fastq = _convert_fastq(fastq_file, pair_file, rg_name,
                                         out_file, config)
            cl = [config["program"].get("mosaik", "MosaikAligner")]
            cl += _mosaik_args_from_config(config)
            cl += extra_args if extra_args is not None else []
            cl += ["-ia", ref_file,
                   "-in", built_fastq,
                   "-out", os.path.splitext(tx_out_file)[0]]
            jump_base = os.path.splitext(ref_file)[0]
            key_file = "{0}_keys.jmp".format(jump_base)
            if file_exists(key_file):
                cl += ["-j", jump_base]
                # XXX hacky way to guess key size which needs to match
                # Can I get hash size directly
                jump_size_gb = os.path.getsize(key_file) / 1073741824.0
                if jump_size_gb < 1.0:
                    cl += ["-hs", "13"]
            cl += _get_mosaik_nn_args(out_file)
            env_set = "export MOSAIK_TMP={0}".format(os.path.dirname(tx_out_file))
            subprocess.check_call(env_set + " && "+
                                  " ".join([str(x) for x in cl]), shell=True)
            os.remove(built_fastq)
    return out_file

def remap_index_fn(ref_file):
    """Map bowtie references to equivalent mosaik indexes.
    """
    return ref_file.replace("/bowtie/", "/mosaik/")
