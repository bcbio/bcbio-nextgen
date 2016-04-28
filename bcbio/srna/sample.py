import os
import sys
import glob
import os.path as op
import shutil
from collections import Counter

try:
    from seqcluster.libs.fastq import collapse, write_output
except ImportError:
    pass

from bcbio.utils import (file_exists, append_stem, replace_directory, symlink_plus)
from bcbio.provenance import do
from bcbio.provenance.versioncheck import java
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio import utils
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.log import logger


def trim_srna_sample(data):
    """
    Remove 3' adapter for smallRNA-seq
    Uses cutadapt but with different parameters than for other pipelines.
    """
    in_file = data["files"][0]
    names = data["rgnames"]['sample']
    work_dir = os.path.join(dd.get_work_dir(data), "trimmed")
    out_dir = os.path.join(work_dir, names)
    utils.safe_makedir(out_dir)
    out_file = replace_directory(append_stem(in_file, ".clean"), out_dir)
    trim_reads = data["config"]["algorithm"].get("trim_reads", True)
    if trim_reads:
        adapter = dd.get_adapters(data)[0]
        out_noadapter_file = replace_directory(append_stem(in_file, ".fragments"), out_dir)
        out_short_file = replace_directory(append_stem(in_file, ".short"), out_dir)
        cutadapt = os.path.join(os.path.dirname(sys.executable), "cutadapt")
        cmd = _cmd_cutadapt()
        if not utils.file_exists(out_file):
            with file_transaction(out_file) as tx_out_file:
                do.run(cmd.format(**locals()), "remove adapter")
    else:
        symlink_plus(in_file, out_file)
    data["clean_fastq"] = out_file
    data["collapse"] = _collapse(data["clean_fastq"])
    data["size_stats"] = _summary(data['collapse'])
    return [[data]]

def sample_annotation(data):
    """
    Annotate miRNAs using miRBase database with seqbuster tool
    """
    names = data["rgnames"]['sample']
    work_dir = os.path.join(dd.get_work_dir(data), "mirbase")
    out_dir = os.path.join(work_dir, names)
    utils.safe_makedir(out_dir)
    out_file = op.join(out_dir, names)
    if dd.get_mirbase_hairpin(data):
        mirbase = op.abspath(op.dirname(dd.get_mirbase_hairpin(data)))
        data['seqbuster'] = _miraligner(data["collapse"], out_file, dd.get_species(data), mirbase, data['config'])
    else:
        logger.debug("No annotation file from miRBase.")

    sps = dd.get_species(data) if dd.get_species(data) else "None"
    if file_exists(op.join(dd.get_work_dir(data), "mirdeep2", "novel", "hairpin.fa")):
        data['seqbuster_novel'] = _miraligner(data["collapse"], "%s_novel" % out_file, sps,  op.join(dd.get_work_dir(data), "mirdeep2", "novel"), data['config'])

    data['trna'] = _trna_annotation(data)
    return [[data]]

def _cmd_cutadapt():
    """
    Run cutadapt for smallRNA data that needs some specific values.
    """
    cmd = "{cutadapt} --adapter={adapter} --untrimmed-output={out_noadapter_file} -o {tx_out_file} -m 17 --overlap=8 {in_file} --too-short-output {out_short_file}"
    return cmd

def _collapse(in_file):
    """
    Collpase reads into unique sequences with seqcluster
    """
    out_file = append_stem(in_file, ".trimming").replace(".gz", "")
    if file_exists(out_file):
        return out_file
    seqs = collapse(in_file)
    write_output(out_file, seqs)
    return out_file

def _summary(in_file):
    """
    Calculate size distribution after adapter removal
    """
    data = Counter()
    out_file = in_file + "_size_stats"
    if file_exists(out_file):
        return out_file
    with open(in_file) as in_handle:
        for line in in_handle:
            counts = int(line.strip().split("_x")[1])
            line = in_handle.next()
            l = len(line.strip())
            in_handle.next()
            in_handle.next()
            data[l] += counts
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, 'w') as out_handle:
            for l, c in data.iteritems():
                out_handle.write("%s %s\n" % (l, c))
    return out_file

def _miraligner(fastq_file, out_file, species, db_folder, config):
    """
    Run miraligner tool (from seqcluster suit) with default
    parameters.
    """
    if not _check_java_version(config, {}):
        return None
    resources = config_utils.get_resources("miraligner", config)
    miraligner = config_utils.get_program("miraligner", config)
    jvm_opts =  "-Xms750m -Xmx4g"
    if resources and resources.get("jvm_opts"):
        jvm_opts = " ".join(resources.get("jvm_opts"))

    cmd = ("{miraligner} {jvm_opts} -freq -sub 1 -trim 3 -add 3 -s {species} -i {fastq_file} -db {db_folder}  -o {tx_out_file}")
    if not file_exists(out_file + ".mirna"):
        with file_transaction(out_file) as tx_out_file:
            do.run(cmd.format(**locals()), "Do miRNA annotation for %s" % fastq_file)
            if _old_version(tx_out_file + ".mirna"):
                raise ValueError("Please install last version for miraligner."
                                 "bcbio_nextgen.py upgrade -u deps --tools.")
            shutil.move(tx_out_file + ".mirna", out_file + ".mirna")
    return out_file + ".mirna"

def _old_version(fn):
    """Check if miraligner is old version."""
    with open(fn) as in_handle:
        h = in_handle.next()
        if h.find("freq") == -1:
            return True
    return False

def _trna_annotation(data):
    """
    use tDRmapper to quantify tRNAs
    """
    trna_ref = op.join(dd.get_srna_trna_file(data))
    name = dd.get_sample_name(data)
    work_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "trna", name))
    in_file = op.basename(data["clean_fastq"])
    tdrmapper = os.path.join(os.path.dirname(sys.executable), "TdrMappingScripts.pl")
    perl_export = utils.get_perl_exports()
    if not file_exists(trna_ref) or not file_exists(tdrmapper):
        logger.info("There is no tRNA annotation to run TdrMapper.")
        return work_dir
    out_file = op.join(work_dir, in_file + ".hq_cs.mapped")
    if not file_exists(out_file):
        with tx_tmpdir(data) as txdir:
            with utils.chdir(txdir):
                utils.symlink_plus(data["clean_fastq"], op.join(txdir, in_file))
                cmd = ("{perl_export} && perl {tdrmapper} {trna_ref} {in_file}").format(**locals())
                do.run(cmd, "tRNA for %s" % name)
                for filename in glob.glob("*mapped*"):
                    shutil.move(filename, work_dir)
    return work_dir

def _check_java_version(config, items):
    msg = java(config, items)
    if msg:
        logger.warning("miraligner is only compatible with java 1.7")
        return False
    return True
