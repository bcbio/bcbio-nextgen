"""Run mirge tool"""
import os
import sys
import shutil
import glob
from bcbio.distributed.transaction import tx_tmpdir
from bcbio.provenance import do
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.log import logger
from bcbio import utils

SPS = {'hsa': 'human',
       'mmu': 'mouse',
       'rno': 'rat',
       'dre': 'zebrafish',
       'cel': 'nematode',
       'dme': 'fruitfly'}

def run(data):
    """Proxy function to run the tool"""
    sample = data[0][0]
    work_dir = dd.get_work_dir(sample)
    out_dir = os.path.join(work_dir, "mirge")
    lib = _find_lib(sample)
    mirge = _find_mirge(sample)
    bowtie = _find_bowtie(sample)
    sps = dd.get_species(sample)
    species = SPS.get(sps, "")
    if not species:
        raise ValueError("species not supported (hsa, mmu, rno, dre, cel, dme): %s" % sps)
    if not lib:
        raise ValueError("-lib option is not set up in resources for mirge tool."
                         " Read above warnings lines.")

    if not utils.file_exists(out_dir):
        with tx_tmpdir() as tmp_dir:
            sample_file = _create_sample_file(data, tmp_dir)
            do.run(_cmd().format(**locals()), "Running miRge2.0.")
            shutil.move(tmp_dir, out_dir)
    return [os.path.abspath(fn) for fn in glob.glob(os.path.join(out_dir, "*", "*"))]

def _find_mirge(data):
    try:
        mirge = config_utils.get_program("miRge2.0", data)
        return mirge
    except config_utils.CmdNotFound:
        logger.warning("miRge2.0 is not found. Install it first, and try again.")
    return None

def _cmd():
    cmd = "{mirge} annotate -s {sample_file} -o {tmp_dir} -pb {bowtie} {lib} -sp {species}   -di -ai -gff"
    return cmd

def _create_sample_file(data, out_dir):
    """from data list all the fastq files in a file"""
    sample_file = os.path.join(out_dir, "sample_file.txt")
    with open(sample_file, 'w') as outh:
        for sample in data:
            outh.write(sample[0]["clean_fastq"] + "\n")
    return sample_file

def _find_bowtie(data):
    try:
        bowtie = config_utils.get_program("bowtie", data)
        return os.path.dirname(bowtie)
    except config_utils.CmdNotFound:
        logger.watning("bowtie is not found. Install it first, and try again.")
    return None

def _find_lib(data):
    """Find mirge libs"""
    options = " ".join(data.get('resources', {}).get('mirge', {}).get("options", ""))
    if options.find("-lib") > -1 and utils.file_exists(options.split()[1]):
        return options
    if not options:
        logger.warning("miRge libraries not found. Follow these instructions to install them:")
        logger.warning("https://github.com/mhalushka/miRge#download-libraries")
        logger.warning("Then, pass -lib LIB_PATH with resourcces:mirge:options:[...]")
        logger.warning("More information: https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#smallrna-seq")
