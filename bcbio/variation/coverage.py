"""Chanjo provides a better way to handle sequence coverage data in
clinical sequencing.  This code integrates it with bcbio
"""

import os.path

import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.utils import safe_makedir

# 
# https://groups.google.com/forum/#!topic/biovalidation/GjFQzJKxl-k 

# <quote> we had a `coverage` input in the `algorithm` section that specified
# a BED file to use for calculating coverage </quote>

# <quote> We could do this or also use the `variant_regions` BED file for 
# assessment. </quote>

def summary (samples, run_parallel):
    cmd_name = 'chanjo'
    
    for data in samples:
        # input 
        bam = tz.get_in(["work_bam"], data[0], None)
        sample_name = tz.get_in(['rgnames','sample'], data[0], None)
        bed_file = tz.get_in(["config", "algorithm", "coverage"], data[0], None)

        outport_dir = os.path.abspath(tz.get_in(['upload', 'dir'], data[0]))
        if not os.path.exists(outport_dir):
            safe_makefir(outport_dir)

        outport = os.path.join(outport_dir, sample_name, '{0}-coverage.bed'.format(sample_name))
        if not utils.file_exists(outport):
            with file_transaction(outport) as tx_out_file:
                cmd = [cmd_name,'annotate', bam, bed_file, "-o", tx_out_file, '-s', sample_name]
                do.run(cmd, "Summarizing coverage with {0}".format(cmd_name), samples[0])
    out = []
    for x in samples[0]:
        outport_dir = os.path.abspath(tz.get_in(['upload', 'dir'], data[0]))
        outport = os.path.join(outport_dir, sample_name, '{0}-coverage.bed'.format(sample_name))
        x["coverage"] = {"summary": outport}
        out.append([x])
    return out

