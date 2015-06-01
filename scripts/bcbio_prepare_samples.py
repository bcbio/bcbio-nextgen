#!/usr/bin/env python -Es
"""
Script that creates bcbio-compatible inputs in case of multiple files samples 
"""


import os
import yaml
from collections import defaultdict
from argparse import ArgumentParser
from cluster_helper import cluster as ipc
from bcbio import log
from bcbio.log import logger
from bcbio.install import _get_data_dir
from bcbio import utils
from bcbio.bam import is_bam
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


def _get_samples_to_process(fn, out_dir, config):
    """parse csv file with one line per file. It will merge
    all files that have the same description name"""
    samples = defaultdict(list)
    with open(fn) as handle:
        for l in handle:
            if not l.startswith("samplename"):
                cols = l.strip().split(",")
                samples[cols[1]].append(cols)
    for sample, items in samples.iteritems():
        if is_fastq(items[0][0], True):
            fn = "fq_merge"
            ext = ".fastq.gz"
        elif is_bam(items[0][0]):
            fn = "bam_merge"
            ext = ".bam"
        files = [os.path.abspath(fn_file[0]) for fn_file in items]
        samples[sample] = [{'files': _check_paired(files), 'out_file': os.path.join(out_dir, sample + ext), 'fn': fn, 'anno': items[0][2:], 'config': config, 'name': sample, 'out_dir': out_dir}]
    return [samples[sample] for sample in samples]


def _check_paired(files):
    """check if files are fastq(.gz) and paired"""
    if files[0].endswith(".bam"):
        return files
    return combine_pairs(files)


def get_cluster_view(p):
    """get ipython running"""
    return ipc.cluster_view(p['scheduler'], p['queue'], p['num_jobs'], p['cores_per_job'], start_wait=p['timeout'], extra_params={"resources": p['resources'], "mem": p['mem'], "tag": p['tag'], "run_local": False})


def wait_until_complete(jobs):
    """wait jobs finish"""
    return [j.get() for j in jobs]


if __name__ == "__main__":
    description = ("Merge multiple files from the same sample to be compatible with bcbio BAM/FASTQ input files")

    parser = ArgumentParser(description="Merge fastq or bam files")
    parser.add_argument("--csv", required=True, help="csv file with metadata")
    parser.add_argument("--out", required=True, help="output dir")
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
    system_config = os.path.join(_get_data_dir(), "galaxy", "bcbio_system.yaml")
    with open(system_config) as in_handle:
        config = yaml.load(in_handle)
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
    samples = _get_samples_to_process(args.csv, out_dir, config)
    parallel = resources.calculate(parallel, [samples], sysinfo, config)

    with prun.start(parallel, samples, config, dirs) as run_parallel:
        with profile.report("prepare bcbio samples", dirs):
            samples = run_parallel("prepare_bcbio_samples", samples)

    create_new_csv(samples, args)
