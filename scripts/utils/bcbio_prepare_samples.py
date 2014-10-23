import os
import yaml
import itertools as tz
from collections import defaultdict
from argparse import ArgumentParser
from cluster_helper import cluster as ipc
from bcbio.log import logger
from bcbio.install import _get_data_dir
from bcbio.bam.fastq import merge as fq_merge
from bcbio.bam import merge as bam_merge
from bcbio.bam import is_bam
from bcbio.bam.fastq import is_fastq
from bcbio.distributed.transaction import file_transaction


def create_new_csv(prep, samples, args):
    """create csv file that can be use with bcbio -w template"""
    ori = os.path.basename(args.csv)
    out_fn = os.path.join(args.out, os.path.splitext(ori)[0] + "-merged.csv")
    logger.info("Preparing new csv: %s" % out_fn)
    with file_transaction(out_fn) as tx_out:
        with open(tx_out, 'w') as handle:
            handle.write(_header(args.csv))
            for p, s in tz.izip(prep, samples.keys()):
                handle.write("%s,%s,%s\n" % (os.path.basename(p), s, ",".join(samples[s]['anno'])))


def _header(fn):
    """read header of csv file"""
    l = open(fn).readline()
    return l


def _get_samples_to_process(fn):
    """parse csv file with one line per file. It will merge
    all files that have the same description name"""
    samples = defaultdict(list)
    with open(fn) as handle:
        for l in handle:
            if not l.startswith("samplename"):
                cols = l.strip().split(",")
                samples[cols[1]].append(cols)
    for sample, anno in samples.iteritems():
        if is_fastq(anno[0][0]):
            fn = fq_merge
            ext = ".fastq.gz"
        elif is_bam(anno[0][0]):
            fn = bam_merge
            ext = ".bam"
        files = [fn_file[0] for fn_file in anno]
        samples[sample] = {'files': files, 'out_file': sample + ext, 'fn': fn, 'anno': anno[0][2:]}
    return samples


def get_cluster_view(args):
    """get ipython running"""
    return ipc.cluster_view(args.scheduler, args.queue, args.num_jobs, args.cores_per_job, start_wait=args.timeout, extra_params={"resources": args.resources,"mem": args.memory_per_job,"tag": "bcbio_prepare","run_local": False})


def wait_until_complete(jobs):
    """wait jobs finish"""
    return [j.get() for j in jobs]


if __name__ == "__main__":
    parser = ArgumentParser(description="Merge fastq or bam files")
    parser.add_argument("--csv", required=True, help="csv file with metadata")
    parser.add_argument("--out", required=True, help="output dir")
    parser.add_argument("-n", "--num-jobs", type=int,
                        default=1, help="Number of concurrent jobs to process.")
    parser.add_argument("-c", "--cores-per-job", type=int,
                        default=1, help="Number of cores to use.")
    parser.add_argument("-m", "--memory-per-job", default=2, help="Memory in GB to reserve per job.")
    parser.add_argument("--timeout", default=15, help="Time to wait before giving up starting.")
    parser.add_argument("-s", "--scheduler", default=None, help="Type of scheduler to use.",
                        choices=["lsf", "slurm", "torque", "sge"])
    parser.add_argument("-r", "--resources", default=None, help="Extra scheduler resource flags.")
    parser.add_argument("-q", "--queue", default=None, help="Queue to submit jobs to.")
    parser.add_argument("-t", "--paralleltype",
                        choices=["local", "ipython"],
                        default="local", help="Run with iptyhon")
    args = parser.parse_args()
    system_config = os.path.join(_get_data_dir(), "galaxy", "bcbio_system.yaml")
    with open(system_config) as in_handle:
        config = yaml.load(in_handle)
        config["algorithm"] = {"num_cores": 1}
    samples = _get_samples_to_process(args.csv)
    prepped = []
    if args.paralleltype == "ipython":
        logger.info("Starting IPython cluster. This may take a while.")
        with get_cluster_view(args) as view:
            logger.info("IPython cluster is up.")
            for sample, info in samples.iteritems():
                prepped.append(view.apply_async(info['fn'], info["files"], os.path.join(args.out, info["out_file"]), config))
            prepped = wait_until_complete(prepped)
    else:
        for sample, info in samples.iteritems():
            logger.info("Merging sample: %s" % sample)
            prepped.append(info['fn'](info["files"], os.path.join(args.out, info["out_file"]), config))
    create_new_csv(prepped, samples, args)
