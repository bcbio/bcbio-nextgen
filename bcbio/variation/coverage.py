"""Chanjo provides a better way to handle sequence coverage data in
clinical sequencing.  This code integrates it with bcbio
"""

import codecs
import os.path

import chanjo
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

def summary(samples, run_parallel):
    cmd_name = 'chanjo'
    for data in samples:
        # input
        bam = tz.get_in(["work_bam"], data[0], None)
        sample_name = tz.get_in(['rgnames', 'sample'], data[0], None)
        bed_file = tz.get_in(["config", "algorithm", "coverage"], data[0], None)

        output_dir = os.path.abspath(tz.get_in(['upload', 'dir'], data[0]))
        if not os.path.exists(output_dir):
            safe_makedir(output_dir)

        output = os.path.join(output_dir, sample_name, '{0}-coverage.bed'.format(sample_name))
        if not utils.file_exists(output):
            with file_transaction(data, output) as tx_out_file:
                with codecs.open(bed_file, encoding='utf-8') as bed_stream:
                    with codecs.open(output, "w", encoding='utf-8') as coverage_stream:
                        for line in chanjo.annotate_bed_stream(bed_stream, bam):
                            coverage_stream.write(chanjo.serialize_interval(line))
                            coverage_stream.write('\n')

    out = []
    for x in samples[0]:
        output_dir = os.path.abspath(tz.get_in(['upload', 'dir'], data[0]))
        output = os.path.join(output_dir, sample_name, '{0}-coverage.bed'.format(sample_name))
        x["coverage"] = {"summary": output}
        out.append([x])
    return out

def summarize_samples(samples, run_parallel):
    """Back compatibility for existing pipelines. Should be replaced with summary when ready.
    """
    return samples
