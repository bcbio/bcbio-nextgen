#!/usr/bin/env python
"""Convert 3 fastq inputs (read 1, read 2, UMI) into paired inputs with UMIs in read names

Usage:
  bcbio_fastq_umi_prep.py single <out basename> <read 1 fastq> <read 2 fastq> <umi fastq>
or:
  bcbio_fastq_umi_prep.py autopair [<list> <of> <fastq> <files>]

Creates two fastq files with embedded UMIs: <out_basename>_R1.fq.gz <out_basename>_R2.fq.gz
or a directory of fastq files with UMIs added to the names.
"""
from __future__ import print_function
import argparse
import os
import sys

from bcbio import utils
from bcbio.bam import fastq
from bcbio.provenance import do
from bcbio.distributed.multi import run_multicore, zeromq_aware_logging

transform_json = r"""{
    "read1": "(?P<name>@.*)\\n(?P<seq>.*)\\n\\+(.*)\\n(?P<qual>.*)\\n",
    "read2": "(?P<name>@.*)\\n(?P<seq>.*)\\n\\+(.*)\\n(?P<qual>.*)\\n",
    "read3": "(@.*)\\n(?P<MB>.*)\\n\\+(.*)\\n(.*)\\n"
}
"""

def run_single(args):
    add_umis_to_fastq(args.out_base, args.read1_fq, args.read2_fq, args.umi_fq, cores=8)

@utils.map_wrap
@zeromq_aware_logging
def add_umis_to_fastq_parallel(out_base, read1_fq, read2_fq, umi_fq, config):
    add_umis_to_fastq(out_base, read1_fq, read2_fq, umi_fq, cores=1)

def add_umis_to_fastq(out_base, read1_fq, read2_fq, umi_fq, cores=1):
    print("Processing", read1_fq, read2_fq, umi_fq)
    out1_fq = out_base + "_R1.fq.gz"
    out2_fq = out_base + "_R2.fq.gz"
    transform_json_file = out_base + "-transform.json"
    with open(transform_json_file, "w") as out_handle:
        out_handle.write(transform_json)
    with utils.open_gzipsafe(read1_fq) as in_handle:
        ex_name = in_handle.readline().split(" ")
        if len(ex_name) == 2:
            fastq_tags_arg = "--keep_fastq_tags"
        else:
            fastq_tags_arg = ""
    cmd = ("umis fastqtransform {fastq_tags_arg} "
           "--fastq1out >(pbgzip -n {cores} -c > {out1_fq}) "
           "--fastq2out >(pbgzip -n {cores} -c > {out2_fq}) "
           "{transform_json_file} {read1_fq} "
           "{read2_fq} {umi_fq}")
    do.run(cmd.format(**locals()), "Add UMIs to paired fastq files")

    os.remove(transform_json_file)

def run_autopair(args):
    outdir = utils.safe_makedir(args.outdir)
    to_run = []
    extras = []
    for fnames in fastq.combine_pairs(sorted(args.files)):
        if len(fnames) == 2:
            to_run.append(fnames)
        elif len(fnames) == 3:
            r1, r2, r3 = sorted(fnames)
            to_run.append([r1, r2])
            extras.append(r3)
        else:
            assert len(fnames) == 1, fnames
            extras.append(fnames[0])
    ready_to_run = []
    for r1, r2 in to_run:
        target = os.path.commonprefix([r1, r2])
        r3 = None
        for test_r3 in extras:
            if (os.path.commonprefix([r1, test_r3]) == target and
                  os.path.commonprefix([r2, test_r3]) == target):
                r3 = test_r3
                break
        assert r3, (r1, r2, extras)
        base_name = os.path.join(outdir, os.path.commonprefix([r1, r2, r3]).rstrip("_R"))
        ready_to_run.append([base_name, r1, r3, r2, {"algorithm": {}, "resources": {}}])

    parallel = {"type": "local", "cores": args.cores, "progs": []}
    run_multicore(add_umis_to_fastq_parallel, ready_to_run, {"algorithm": {}}, parallel)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add UMIs to fastq read names")
    sp = parser.add_subparsers(title="[sub-commands]")

    p = sp.add_parser("autopair", help="Automatically pair R1/R2/R3 fastq inputs")
    p.add_argument("-c", "--cores", default=1, type=int, help="Number of cores to run in parallel")
    p.add_argument("--outdir", default="with_umis", help="Output directory to write UMI prepped fastqs")
    p.add_argument("files", nargs="*", help="All fastq files to pair and process")
    p.set_defaults(func=run_autopair)

    p = sp.add_parser("single", help="Run single set of fastq files with separate UMIs")
    p.add_argument("out_base", help="Base name for output files -- you get <base_name>_R1.fq.gz")
    p.add_argument("read1_fq", help="Input fastq, read 1")
    p.add_argument("read2_fq", help="Input fastq, read 2")
    p.add_argument("umi_fq", help="Input fastq, UMIs")
    p.set_defaults(func=run_single)
    if len(sys.argv) == 1:
        parser.print_help()
    args = parser.parse_args()
    args.func(args)
