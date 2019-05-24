#!/usr/bin/env python
"""Convert fastq inputs into paired inputs with UMIs in read names.

Handles two cases:

- Separate UMI read files (read 1, read 2, UMI)

  Usage:
    bcbio_fastq_umi_prep.py single <out basename> <read 1 fastq> <read 2 fastq> <umi fastq>
  or:
    bcbio_fastq_umi_prep.py autopair [<list> <of> <fastq> <files>]

- Duplex barcodes where the tags are incorporated into read 1 and read 2

  Usage:
    bcbio_fastq_umi_prep.py single --tag1 5 --tag2 5 <out basename> <read 1 fastq> <read 2 fastq>
  or:
    bcbio_fastq_umi_prep.py autopair --tag1 5 --tag2 5 [<list> <of> <fastq> <files>]

Creates two fastq files with embedded UMIs: <out_basename>_R1.fq.gz <out_basename>_R2.fq.gz
or a directory of fastq files with UMIs added to the names.

autopair assumes 3 sets of reads based on the way bcl2fastq exports index reads:

- R1 -- The first read pair
- R2 -- The UMI barcode reads
- R3 -- The second read pair

If you're using a different approach to generate the UMIs, please
maintain the same R1/R2/R3 naming scheme.
"""
from __future__ import print_function
import argparse
import math
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

duplex_transform = r"""{
    "read1": "(?P<name>@.*)\\n(?P<CB1>.{%s})(?P<seq>.*)\\n\\+(.*)\\n(.{%s})(?P<qual>.*)\\n",
    "read2": "(?P<name>@.*)\\n(?P<CB2>.{%s})(?P<seq>.*)\\n\\+(.*)\\n(.{%s})(?P<qual>.*)\\n"
}"""

r1_transform = r"""{
    "read1": "(?P<name>@.*)\\n(?P<MB>.{%s})(?P<seq>.*)\\n\\+(.*)\\n(.{%s})(?P<qual>.*)\\n",
    "read2": "(?P<name>@.*)\\n(?P<seq>.*)\\n\\+(.*)\\n(?P<qual>.*)\\n",
}"""

r2_transform = r"""{
    "read1": "(?P<name>@.*)\\n(?P<seq>.*)\\n\\+(.*)\\n(?P<qual>.*)\\n",
    "read2": "(?P<name>@.*)\\n(?P<MB>.{%s})(?P<seq>.*)\\n\\+(.*)\\n(.{%s})(?P<qual>.*)\\n"
}"""

def run_single(args):
    tags = [args.tag1, args.tag2] if args.tag1 and args.tag2 else None
    if args.umi_fq and not tags:
        add_umis_to_fastq(args.out_base, args.read1_fq, None, args.umi_fq, tags, cores=args.cores)
        add_umis_to_fastq(args.out_base, None, args.read2_fq, args.umi_fq, tags, cores=args.cores)
    else:
        add_umis_to_fastq(args.out_base, args.read1_fq, args.read2_fq, args.umi_fq, tags, cores=args.cores)

@utils.map_wrap
@zeromq_aware_logging
def add_umis_to_fastq_parallel(out_base, read1_fq, read2_fq, umi_fq, tags, cores, config):
    add_umis_to_fastq(out_base, read1_fq, read2_fq, umi_fq, tags, cores)

def add_umis_to_fastq(out_base, read1_fq, read2_fq, umi_fq, tags=None, cores=1):
    print("Adding UMIs from", umi_fq, "to read headers in", " ".join([x for x in [read1_fq, read2_fq] if x]))
    out1_fq = out_base + "_R1.fq.gz"
    out2_fq = out_base + "_R2.fq.gz"
    if umi_fq and not tags:
        if read1_fq:
            assert not read2_fq
            return _add_umis_with_fastp(read1_fq, umi_fq, out1_fq, cores)
        else:
            assert read2_fq
            return _add_umis_with_fastp(read2_fq, umi_fq, out2_fq, cores)

    transform_json_file = out_base + "-transform.json"
    with open(transform_json_file, "w") as out_handle:
        if tags:
            tag1, tag2 = tags
            if tag1 and tag2:
                out_handle.write(duplex_transform % (tag1, tag1, tag2, tag2))
            elif tag1:
                out_handle.write(r1_transform % (tag1, tag1))
            else:
                out_handle.write(r2_transform % (tag2, tag2))
        else:
            out_handle.write(transform_json)
    with utils.open_gzipsafe(read1_fq) as in_handle:
        ex_name = in_handle.readline().split(" ")
        fastq_tags_arg = "--keep_fastq_tags" if len(ex_name) == 2 else ""
    tag_arg = "--separate_cb" if tags else ""
    cmd = ("umis fastqtransform {fastq_tags_arg} {tag_arg} "
           "--fastq1out >(bgzip --threads {cores} -c > {out1_fq}) "
           "--fastq2out >(bgzip --threads {cores} -c > {out2_fq}) "
           "{transform_json_file} {read1_fq} "
           "{read2_fq}")
    if umi_fq:
        cmd += " {umi_fq}"
    do.run(cmd.format(**locals()), "Add UMIs to paired fastq files")

    os.remove(transform_json_file)

def _add_umis_with_fastp(read_fq, umi_fq, out_fq, cores):
    """Add UMIs to reads from separate UMI file using fastp.
    """
    with utils.open_gzipsafe(umi_fq) as in_handle:
        in_handle.readline()  # name
        umi_size = len(in_handle.readline().strip())
    cmd = ("fastp -Q -A -L -G -w 1 --in1 {read_fq} --in2 {umi_fq} "
           "--umi --umi_prefix UMI --umi_loc read2 --umi_len {umi_size} "
           "--out1 >(bgzip --threads {cores} -c > {out_fq}) --out2 /dev/null "
           "-j /dev/null -h /dev/null")
    do.run(cmd.format(**locals()), "Add UMIs to fastq file with fastp")

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
    tags = [args.tag1, args.tag2] if args.tag1 and args.tag2 else None
    if not tags:
        # Aim for 2 or 3 simultaneous processes, each with multiple cores
        target_processes = 2
        process_cores = max(1, (args.cores // target_processes) + (args.cores % target_processes))
        overall_processes = max(1, int(math.ceil(args.cores / float(process_cores))))
    else:
        process_cores = 1
        overall_processes = args.cores
    for r1, r2 in to_run:
        target = _commonprefix([r1, r2])
        if tags:
            base_name = os.path.join(outdir, os.path.basename(_commonprefix([r1, r2])))
            umi = None
        else:
            r3 = None
            for test_r3 in extras:
                if (_commonprefix([r1, test_r3]) == target and _commonprefix([r2, test_r3]) == target):
                    r3 = test_r3
                    break
            assert r3, (r1, r2, extras)
            base_name = os.path.join(outdir, os.path.basename(_commonprefix([r1, r2, r3])))
            r1, r2, umi = _find_umi([r1, r2, r3])
        # fastp handles a single pair of reads so we split processing to run on each
        if umi and not tags:
            ready_to_run.append([base_name, r1, None, umi, None, process_cores, {"algorithm": {}, "resources": {}}])
            ready_to_run.append([base_name, None, r2, umi, None, process_cores, {"algorithm": {}, "resources": {}}])
        else:
            ready_to_run.append([base_name, r1, r2, umi, tags, process_cores, {"algorithm": {}, "resources": {}}])
    parallel = {"type": "local", "cores": overall_processes, "progs": []}
    run_multicore(add_umis_to_fastq_parallel, ready_to_run, {"algorithm": {}}, parallel)

def _find_umi(files):
    """Find UMI file using different naming schemes.

    R1/R2/R3 => R1/R3 with R2 UMI
    R1/R2/I1 => R1/R2 with I1 UMI
    """
    base = os.path.basename(_commonprefix(files))

    def _file_ext(f):
        exts = utils.splitext_plus(os.path.basename(f).replace(base, ""))[0].split("_")
        exts = [x for x in exts if x]
        return exts[0]

    exts = dict([(_file_ext(f), f) for f in files])
    if "I1" in exts:
        return exts["R1"], exts["R2"], exts["I1"]
    else:
        assert "R3" in exts, exts
        return exts["R1"], exts["R3"], exts["R2"]

def _commonprefix(files):
    """Retrieve a common prefix for files without extra _R1 _I1 extensions.

    Allows alternative naming schemes (R1/R2/R3) (R1/R2/I1).
    """
    out = os.path.commonprefix(files)
    out = out.rstrip("_R")
    out = out.rstrip("_I")
    out = out.rstrip("_")
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add UMIs to fastq read names")
    sp = parser.add_subparsers(title="[sub-commands]")

    p = sp.add_parser("autopair", help="Automatically pair R1/R2 (and maybe R3) fastq inputs")
    p.add_argument("-c", "--cores", default=1, type=int, help="Number of cores, allowing running samples in parallel")
    p.add_argument("--outdir", default="with_umis", help="Output directory to write UMI prepped fastqs")
    p.add_argument("--tag1", help="Duplex read 1 tag -- bases to trim from 5' end")
    p.add_argument("--tag2", help="Duplex read 2 tag -- bases to trim from 5' end")
    p.add_argument("files", nargs="*", help="All fastq files to pair and process")
    p.set_defaults(func=run_autopair)

    p = sp.add_parser("single", help="Run single set of fastq files with UMIs/duplexes")
    p.add_argument("-c", "--cores", default=1, type=int, help="Number of cores to use for parallel bgzip")
    p.add_argument("--tag1", help="Duplex read 1 tag -- bases to trim from 5' end", type=int)
    p.add_argument("--tag2", help="Duplex read 2 tag -- bases to trim from 5' end", type=int)
    p.add_argument("out_base", help="Base name for output files -- you get <base_name>_R1.fq.gz")
    p.add_argument("read1_fq", help="Input fastq, read 1")
    p.add_argument("read2_fq", help="Input fastq, read 2")
    p.add_argument("umi_fq", help="Input fastq, UMIs", nargs="?")
    p.set_defaults(func=run_single)
    if len(sys.argv) == 1:
        parser.print_help()
    args = parser.parse_args()
    args.func(args)
