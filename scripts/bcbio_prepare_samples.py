#!/usr/bin/env python -Es
"""
Script that creates bcbio-compatible inputs in case of multiple files samples
"""


import os
import sys
import yaml
from collections import defaultdict
from argparse import ArgumentParser
from bcbio import log
from bcbio.log import logger
from bcbio.install import _get_data_dir
from bcbio import utils
from bcbio.bam import is_bam
from bcbio.pipeline.sra import is_gsm, is_srr
from bcbio.bam.fastq import is_fastq, combine_pairs
from bcbio.distributed.transaction import file_transaction
from bcbio.distributed import clargs, resources, prun
from bcbio.provenance import system, profile

def create_new_csv(samples, args):
    """create csv file that can be use with bcbio -w template"""
    out_fn = os.path.splitext(args.csv)[0] + "-merged.csv"
    logger.info("Preparing new csv: %s" % out_fn)
    with file_transaction(out_fn) as tx_out:
        with open(tx_out, 'w') as handle:
            handle.write(_header(args.csv))
            for s in samples:
                sample_name = s['name'] if isinstance(s['out_file'], list) else os.path.basename(s['out_file'])
                handle.write("%s,%s,%s\n" % (sample_name, s['name'], ",".join(s['anno'])))


def _header(fn):
    """read header of csv file"""
    l = open(fn).readline()
    return l


def _get_samples_to_process(fn, out_dir, config, force_single, separators):
    """parse csv file with one line per file. It will merge
    all files that have the same description name"""
    out_dir = os.path.abspath(out_dir)
    samples = defaultdict(list)
    with open(fn) as handle:
        for l in handle:
            if l.find("description") > 0:
                logger.info("Skipping header.")
                continue
            cols = l.strip().split(",")
            if len(cols) > 0:
                if len(cols) < 2:
                    raise ValueError("Line needs 2 values: file and name.")
                if utils.file_exists(cols[0]) or is_gsm(cols[0]) or is_srr(cols[0]):
                    if cols[0].find(" ") > -1:
                        new_name = os.path.abspath(cols[0].replace(" ", "_"))
                        logger.warning("Space finds in %s. Linked to %s." % (cols[0], new_name))
                        logger.warning("Please, avoid names with spaces in the future.")
                        utils.symlink_plus(os.path.abspath(cols[0]), new_name)
                        cols[0] = new_name
                    samples[cols[1]].append(cols)
                else:
                    logger.info("skipping %s, File doesn't exist." % cols[0])
    for sample, items in samples.items():
        if is_fastq(items[0][0], True):
            fn = "fq_merge"
            ext = ".fastq.gz"
        elif is_bam(items[0][0]):
            fn = "bam_merge"
            ext = ".bam"
        elif is_gsm(items[0][0]):
            fn = "query_gsm"
            ext = ".fastq.gz"
        elif is_srr(items[0][0]):
            fn = "query_srr"
            ext = ".fastq.gz"
        files = [os.path.abspath(fn_file[0]) if utils.file_exists(fn_file[0]) else fn_file[0] for fn_file in items]
        samples[sample] = [{'files': _check_paired(files, force_single, separators),
                            'out_file': os.path.join(out_dir, sample + ext),
                            'fn': fn, 'anno': items[0][2:], 'config': config,
                            'name': sample, 'out_dir': out_dir}]
    return [samples[sample] for sample in samples]


def _check_stems(files):
    """check if stem names are the same and use full path then"""
    used = set()
    for fn in files:
        if os.path.basename(fn) in used:
            logger.warning("%s stem is multiple times in your file list, "
                         "so we don't know "
                         "how to assign it to the sample data in the CSV. "
                         "We are gonna use full path to make a difference, "
                         "that means paired files should be in the same folder. "
                         "If this is a problem, you should rename the files you want "
                         "to merge. Sorry, no possible magic here." % os.path.basename(fn)
                         )
            return True
        used.add(os.path.basename(fn))
    return False


def _check_paired(files, force_single, separators):
    """check if files are fastq(.gz) and paired"""
    full_name = _check_stems(files)
    if files[0].endswith(".bam"):
        return files
    elif is_gsm(files[0]):
        return files
    return combine_pairs(files, force_single, full_name, separators)


def get_cluster_view(p):
    """get ipython running"""
    from cluster_helper import cluster as ipc
    return ipc.cluster_view(p['scheduler'], p['queue'], p['num_jobs'], p['cores_per_job'], start_wait=p['timeout'], extra_params={"resources": p['resources'], "mem": p['mem'], "tag": p['tag'], "run_local": False})


def wait_until_complete(jobs):
    """wait jobs finish"""
    return [j.get() for j in jobs]


if __name__ == "__main__":
    description = ("Merge multiple files from the same sample to be compatible with bcbio BAM/FASTQ input files")

    parser = ArgumentParser(description="Merge fastq or bam files")
    parser.add_argument("--csv", required=True, help="csv file with metadata")
    parser.add_argument("--out", required=True, help="output dir")
    parser.add_argument("--force-single", action='store_true', default=False, help="Treat all files as single reads")
    parser.add_argument("--separators", nargs="*",
                        default=["R", "_", "-", "."],
                        help="Space separated list of separators that indicates paired files.")
    parser.add_argument("--remove-source", action='store_true', default=False,
                        help="Remove original files.")
    parser.add_argument("-n", "--numcores", type=int,
                        default=1, help="Number of concurrent jobs to process.")
    parser.add_argument("-c", "--cores-per-job", type=int,
                        default=1, help="Number of cores to use.")
    parser.add_argument("-m", "--memory-per-job", default=2, help="Memory in GB to reserve per job.")
    parser.add_argument("--timeout", default=15, help="Time to wait before giving up starting.")
    parser.add_argument("--retries", default=0, type=int,
                        help=("Number of retries of failed tasks during "
                              "distributed processing. Default 0 "
                              "(no retries)"))
    parser.add_argument("-s", "--scheduler", help="Type of scheduler to use.",
                        choices=["lsf", "slurm", "torque", "sge", "pbspro"])
    parser.add_argument("-r", "--resources", help="Extra scheduler resource flags.", default=[], action="append")
    parser.add_argument("-q", "--queue", help="Queue to submit jobs to.")
    parser.add_argument("-p", "--tag", help="Tag name to label jobs on the cluster", default="bcb-prep")
    parser.add_argument("-t", "--paralleltype",
                        choices=["local", "ipython"],
                        default="local", help="Run with iptyhon")

    args = parser.parse_args()
    out_dir = os.path.abspath(args.out)
    utils.safe_makedir(out_dir)
    try:
        system_config = os.path.join(_get_data_dir(), "galaxy", "bcbio_system.yaml")
    except ValueError as err:
        print(err)
        print("WARNING: Attempting to read bcbio_system.yaml in the current directory.")
        system_config = "bcbio_system.yaml"

    if utils.file_exists(system_config):
        with open(system_config) as in_handle:
            config = yaml.safe_load(in_handle)
    else:
        print("WARNING: bcbio_system.yaml not found, creating own resources.")
        config = {'resources': {}}
    res = {'cores': args.cores_per_job}
    config["algorithm"] = {"num_cores": args.cores_per_job}
    config["resources"].update({'sambamba': res,
                                'samtools': res})
    config["log_dir"] = os.path.join(os.path.abspath(os.getcwd()), "log")

    parallel = clargs.to_parallel(args)
    parallel.update({'progs': ['samtools', 'sambamba']})
    parallel = log.create_base_logger(config, parallel)
    log.setup_local_logging(config, parallel)
    dirs = {'work': os.path.abspath(os.getcwd())}
    system.write_info(dirs, parallel, config)
    sysinfo = system.machine_info()[0]
    config["remove_source"] = args.remove_source
    samples = _get_samples_to_process(args.csv, out_dir, config, args.force_single, args.separators)
    if not samples:
        print("No samples found.")
        sys.exit(0)
    parallel = resources.calculate(parallel, [samples], sysinfo, config)

    with prun.start(parallel, samples, config, dirs) as run_parallel:
        with profile.report("prepare bcbio samples", dirs):
            samples = run_parallel("prepare_bcbio_samples", samples)

    create_new_csv(samples, args)
