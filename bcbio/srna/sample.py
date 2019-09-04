import os
import sys
import glob
import re
import os.path as op
import shutil
from collections import Counter

try:
    from dnapilib.apred import iterative_adapter_prediction
    error_dnapi = None
except ImportError:
    error_dnapi = ("No dnapi installed. Need to give adapter sequence."
                   "Please, install with bcbio_conda install dnapi -c bioconda"
                   " or add adapters: ['ADAPTER_SEQ'] to config file.")
    pass

from bcbio.utils import (file_exists, append_stem, replace_directory, symlink_plus)
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio import utils
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.log import logger
from bcbio.rnaseq import spikein
from bcbio.srna.umis import umi_transform


def trim_srna_sample(data):
    """
    Remove 3' adapter for smallRNA-seq
    Uses cutadapt but with different parameters than for other pipelines.
    """
    data = umi_transform(data)
    in_file = data["files"][0]
    names = data["rgnames"]['sample']
    work_dir = os.path.join(dd.get_work_dir(data), "trimmed")
    out_dir = os.path.join(work_dir, names)
    log_out = os.path.join(out_dir, "%s.log" % names)
    utils.safe_makedir(out_dir)
    out_file = replace_directory(append_stem(in_file, ".clean"), out_dir)
    trim_reads = data["config"]["algorithm"].get("trim_reads", True)
    if utils.file_exists(out_file):
        data["files"][0] = out_file
        data["clean_fastq"] = out_file
        data["collapse"] = _collapse(data["clean_fastq"])
        data["size_stats"] = _summary(data['collapse'])
        data["log_trimming"] = log_out
        return [[data]]

    adapter = dd.get_adapters(data)
    is_4n = any([a == "4N" for a in adapter])
    adapter = [a for a in adapter if re.compile(r"^([NATGC]+)$").match(a)]
    if adapter and not trim_reads:
        trim_reads = True
        logger.info("Adapter is set up in config file, but trim_reads is not true."
                    "If you want to skip trimming, skip adapter option from config.")
    if trim_reads and not adapter and error_dnapi:
        raise ValueError(error_dnapi)
    if trim_reads:
        adapters = adapter if adapter else _dnapi_prediction(in_file, out_dir)
    times = "" if not trim_reads or len(adapters) == 1 else "--times %s" % len(adapters)
    if trim_reads and adapters:
        adapter_cmd = " ".join(map(lambda x: "-a " + x, adapters))
        if any([a for a in adapters if re.compile(r"^N+$").match(a)]):
            adapter_cmd = "-N %s" % adapter_cmd
        out_noadapter_file = replace_directory(append_stem(in_file, ".fragments"), out_dir)
        out_short_file = replace_directory(append_stem(in_file, ".short"), out_dir)
        # atropos = _get_atropos()
        atropos = config_utils.get_program("atropos", data, default="atropos")
        options = " ".join(data.get('resources', {}).get('atropos', {}).get("options", ""))
        if options.strip() == "-u 4 -u -4":
            options = ""
            is_4n = "4N"
        cores = ("--threads %s" % dd.get_num_cores(data) if dd.get_num_cores(data) > 1 else "")
        if " ".join(data.get('resources', {}).get('cutadapt', {}).get("options", "")):
            raise ValueError("Atropos is now used, but cutadapt options found in YAML file."
                             "See https://atropos.readthedocs.io/en/latest/")
        cmd = _cmd_atropos()
        if not utils.file_exists(out_file):
            with file_transaction(out_file) as tx_out_file:
                do.run(cmd.format(**locals()), "remove adapter for %s" % names)
                if utils.file_exists(log_out):
                    content = open(log_out).read().replace(out_short_file, names)
                    open(log_out, 'w').write(content)
                if is_4n:
                    options = "-u 4 -u -4"
                    in_file = append_stem(tx_out_file, ".tmp")
                    utils.move_safe(tx_out_file, in_file)
                    cmd = "{atropos} {cores} {options} -se {in_file} -o {tx_out_file} -m 17"
                    do.run(cmd.format(**locals()), "atropos with this parameters %s for %s" %(options, names))
        data["log_trimming"] = log_out
    else:
        if not trim_reads:
            logger.debug("Skip trimming for: %s" % names)
        elif not adapters:
            logger.info("No adapter founds in %s, this is an issue related"
                        " to no small RNA enrichment in your sample." % names)
        symlink_plus(in_file, out_file)
    data["files"][0] = out_file
    data["clean_fastq"] = out_file
    data["collapse"] = _collapse(data["clean_fastq"])
    data["size_stats"] = _summary(data['collapse'])
    return [[data]]

def sample_annotation(data):
    """
    Annotate miRNAs using miRBase database with seqbuster tool
    """
    names = data["rgnames"]['sample']
    tools = dd.get_expression_caller(data)
    work_dir = os.path.join(dd.get_work_dir(data), "mirbase")
    out_dir = os.path.join(work_dir, names)
    utils.safe_makedir(out_dir)
    out_file = op.join(out_dir, names)
    if dd.get_mirbase_hairpin(data):
        mirbase = op.abspath(op.dirname(dd.get_mirbase_hairpin(data)))
        if utils.file_exists(data["collapse"]):
            data['transcriptome_bam'] = _align(data["collapse"],
                                               dd.get_mirbase_hairpin(data),
                                               out_file,
                                               data)
            data['seqbuster'] = _miraligner(data["collapse"], out_file,
                                            dd.get_species(data),
                                            mirbase,
                                            data['config'])
            data["mirtop"] = _mirtop(data['seqbuster'],
                                     dd.get_species(data),
                                     mirbase,
                                     out_dir,
                                     data['config'])
        else:
            logger.debug("Trimmed collapsed file is empty for %s." % names)
    else:
        logger.debug("No annotation file from miRBase.")

    sps = dd.get_species(data) if dd.get_species(data) else "None"
    logger.debug("Looking for mirdeep2 database for %s" % names)
    if file_exists(op.join(dd.get_work_dir(data), "mirdeep2", "novel", "hairpin.fa")):
        data['seqbuster_novel'] = _miraligner(data["collapse"], "%s_novel" % out_file, sps,
                                              op.join(dd.get_work_dir(data),
                                                      "mirdeep2", "novel"),
                                              data['config'])

    if "trna" in tools:
        data['trna'] = _mint_trna_annotation(data)

    data = spikein.counts_spikein(data)
    return [[data]]

def _prepare_file(fn, out_dir):
    """Cut the beginning of the reads to avoid detection of miRNAs"""
    atropos = _get_atropos()
    cmd = "{atropos} trim --max-reads 500000 -u 22 -se {fn} -o {tx_file}"
    out_file = os.path.join(out_dir, append_stem(os.path.basename(fn), "end"))
    if file_exists(out_file):
        return out_file
    with file_transaction(out_file) as tx_file:
        do.run(cmd.format(**locals()))
    return out_file

def _dnapi_prediction(fn, out_dir):
    end_file = _prepare_file(fn, out_dir)
    iterative_result = iterative_adapter_prediction(end_file, [1.2, 1.3, 1.4, 1.7, 2], [7, 11], 500000)
    max_score = iterative_result[1][1]
    adapters = list()
    for a in iterative_result:
        if a[1] > max_score * 0.40:
            logger.debug("Adding adapter to the list: %s with score %s" % (a[0], a[1]))
            adapters.append(a[0])
    return adapters

def _cmd_atropos():
    """
    Run cutadapt for smallRNA data that needs some specific values.
    """
    cmd = "{atropos} {cores} {times} {options} {adapter_cmd} --untrimmed-output={out_noadapter_file} -o {tx_out_file} -m 17 --overlap=8 -se {in_file} --too-short-output {out_short_file} | tee > {log_out}"
    return cmd

def _collapse(in_file):
    """
    Collpase reads into unique sequences with seqcluster
    """
    seqcluster = op.join(utils.get_bcbio_bin(), "seqcluster")
    out_file = "%s.fastq" % utils.splitext_plus(append_stem(in_file, "_trimmed"))[0]
    out_dir = os.path.dirname(in_file)
    if file_exists(out_file):
        return out_file
    cmd = ("{seqcluster} collapse -o {out_dir} -f {in_file} -m 1 --min_size 16")
    do.run(cmd.format(**locals()), "Running seqcluster collapse in %s." % in_file)
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
            line = next(in_handle)
            l = len(line.strip())
            next(in_handle)
            next(in_handle)
            data[l] += counts
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, 'w') as out_handle:
            for l, c in data.items():
                out_handle.write("%s %s\n" % (l, c))
    return out_file

def _align(infile, ref, out_file, data):
    out_file = "%s.bam" % out_file
    razers3 = config_utils.get_program("razers3", data["config"])
    cmd = "{razers3} -dr 0 -i 80 -rr 90 -f -o {tx_out} {ref} {infile}"
    if not razers3:
        logger.info("razers3 is not installed, skipping BAM file creation")
        return None
    if not file_exists(out_file):
        with file_transaction(data, out_file) as tx_out:
            do.run(cmd.format(**locals()), "Running razers3 against hairpins with %s" % infile)
    return out_file

def _miraligner(fastq_file, out_file, species, db_folder, config):
    """
    Run miraligner tool (from seqcluster suit) with default
    parameters.
    """
    resources = config_utils.get_resources("miraligner", config)
    miraligner = config_utils.get_program("miraligner", config)
    jvm_opts =  "-Xms750m -Xmx4g"
    if resources and resources.get("jvm_opts"):
        jvm_opts = " ".join(resources.get("jvm_opts"))
    export = _get_env()
    cmd = ("{export} {miraligner} {jvm_opts} -freq -sub 1 -trim 3 -add 3 -minl 16"
           " -s {species} -i {fastq_file} -db {db_folder}  -o {tx_out_file}")
    if not file_exists(out_file + ".mirna"):
        with file_transaction(out_file) as tx_out_file:
            do.run(cmd.format(**locals()), "Do miRNA annotation for %s" % fastq_file)
            shutil.move(tx_out_file + ".mirna", out_file + ".mirna")
    return out_file + ".mirna"

def _get_env():
    anaconda_bin = os.path.dirname(utils.Rscript_cmd())
    return "unset JAVA_HOME  && export PATH=%s:\"$PATH\" && " % (anaconda_bin)

def _get_atropos():
    anaconda = os.path.dirname(os.path.realpath(sys.executable))
    return os.path.join(anaconda, "..", "envs", "python3", "bin", "atropos")

def _mirtop(input_fn, sps, db, out_dir, config):
    """
    Convert to GFF3 standard format
    """
    hairpin = os.path.join(db, "hairpin.fa")
    gtf = os.path.join(db, "mirbase.gff3")
    if not file_exists(hairpin) or not file_exists(gtf):
        logger.warning("%s or %s are not installed. Skipping." % (hairpin, gtf))
        return None
    out_gtf_fn = "%s.gtf" % utils.splitext_plus(os.path.basename(input_fn))[0]
    out_gff_fn = "%s.gff" % utils.splitext_plus(os.path.basename(input_fn))[0]
    export = _get_env()
    cmd = ("{export} mirtop gff  --sps {sps} --hairpin {hairpin} "
           "--gtf {gtf} --format seqbuster -o {out_tx} {input_fn}")
    if not file_exists(os.path.join(out_dir, out_gtf_fn)) and \
       not file_exists(os.path.join(out_dir, out_gff_fn)):
        with tx_tmpdir() as out_tx:
            do.run(cmd.format(**locals()), "Do miRNA annotation for %s" % input_fn)
            with utils.chdir(out_tx):
                out_fn = out_gtf_fn if utils.file_exists(out_gtf_fn) \
                                    else out_gff_fn
                if utils.file_exists(out_fn):
                    shutil.move(os.path.join(out_tx, out_fn),
                                os.path.join(out_dir, out_fn))
    out_fn = out_gtf_fn if utils.file_exists(os.path.join(out_dir, out_gtf_fn)) \
                        else os.path.join(out_dir, out_gff_fn)
    if utils.file_exists(os.path.join(out_dir, out_fn)):
        return os.path.join(out_dir, out_fn)

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

def _mint_trna_annotation(data):
    """
    use MINTmap to quantify tRNAs
    """
    name = dd.get_sample_name(data)
    work_dir = os.path.join(dd.get_work_dir(data), "trna_mint", name)
    if not dd.get_srna_mint_lookup(data):
        logger.info("There is no tRNA annotation to run MINTmap.")
        return work_dir
    trna_lookup = op.join(dd.get_srna_mint_lookup(data))
    trna_space = op.join(dd.get_srna_mint_space(data))
    trna_other = op.join(dd.get_srna_mint_other(data))
    in_file = op.basename(data["clean_fastq"])
    mintmap = os.path.realpath(os.path.join(os.path.dirname(sys.executable), "MINTmap.pl"))
    perl_export = utils.get_perl_exports()
    if not file_exists(trna_lookup) or not file_exists(mintmap):
        logger.info("There is no tRNA annotation to run MINTmap.")
        return work_dir
    jar_folder = os.path.join(os.path.dirname(mintmap), "MINTplates")
    out_file = op.join(work_dir, name + "-MINTmap_v1-exclusive-tRFs.expression.txt")
    if not file_exists(out_file):
        with tx_tmpdir(data) as txdir:
            with utils.chdir(txdir):
                utils.symlink_plus(data["clean_fastq"], op.join(txdir, in_file))
                cmd = ("{perl_export} && {mintmap} -f {in_file} -p {name} "
                       "-l {trna_lookup} -s {trna_space} -j {jar_folder} "
                       "-o {trna_other}").format(**locals())
                do.run(cmd, "tRNA for %s" % name)
                for filename in glob.glob("*MINTmap*"):
                    shutil.move(filename, work_dir)
    return work_dir
