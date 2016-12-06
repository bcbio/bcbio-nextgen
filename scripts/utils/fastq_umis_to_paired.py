#!/usr/bin/env python
"""Convert 3 fastq inputs (read 1, read 2, UMI) into paired inputs with UMIs in read names

Usage:
  fastq_umis_to_paired.py <out basename> <read 1 fastq> <read 2 fastq> <umi fastq>

Creates two fastq files with embedded UMIs: <out_basename>_R1.fq.gz <out_basename>_R2.fq.gz
"""

transform_json = r"""{
    "read1": "(?P<name>@.*)\\n(?P<seq>.*)\\n\\+(.*)\\n(?P<qual>.*)\\n",
    "read2": "(?P<name>@.*)\\n(?P<seq>.*)\\n\\+(.*)\\n(?P<qual>.*)\\n",
    "read3": "(@.*)\\n(?P<MB>.*)\\n\\+(.*)\\n(.*)\\n"
}
"""
import os
import sys

from bcbio import utils
from bcbio.provenance import do

def main(out_base, read1_fq, read2_fq, umi_fq):
    cores = 8
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
    cmd = ("umis fastqtransform {fastq_tags_arg} --umi_only "
           "--fastq1out >(pbgzip -n {cores} -c > {out1_fq}) "
           "--fastq2out >(pbgzip -n {cores} -c > {out2_fq}) "
           "{transform_json_file} {read1_fq} "
           "{read2_fq} {umi_fq}")
    do.run(cmd.format(**locals()), "Add UMIs to paired fastq files")

    os.remove(transform_json_file)

if __name__ == "__main__":
    main(*sys.argv[1:])
