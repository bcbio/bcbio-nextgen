"""
Unique Molecular Identifier (UMI) handling for small RNAseq.
Most of this either uses Valentine Svennson's umis repository or adapts
code written from it.
https://github.com/vals/umis
"""
import os
import sys
import glob
import re

import bcbio.pipeline.datadict as dd
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.utils import (file_exists, safe_makedir)
from bcbio.distributed.transaction import file_transaction
from bcbio.bam.fastq import open_fastq
from bcbio.rnaseq.umi import get_transform_file, \
                             is_transformed
from bcbio.log import logger

TD = os.path.join(os.path.dirname(__file__), os.pardir, "data", "umis")
p = re.compile(r'smallrna', re.IGNORECASE)

TF = glob.glob(os.path.join(TD, "*-transform.json"))
TF = [t for t in TF if p.match(t)]

SUPPORTED_TRANSFORMS = set([os.path.basename(x).replace("-transform.json", "")
                            for x in TF])


def umi_transform(data):
    """
    transform each read by identifying the barcode and UMI for each read
    and putting the information in the read name
    """
    fq1 = data["files"][0]
    umi_dir = os.path.join(dd.get_work_dir(data), "umis")
    safe_makedir(umi_dir)
    transform = dd.get_umi_type(data)

    if not transform:
        logger.info("No UMI transform specified, assuming pre-transformed data.")
        if is_transformed(fq1):
            logger.info("%s detected as pre-transformed, passing it on unchanged." % fq1)
            data["files"] = [fq1]
            return data
        else:
            logger.error("No UMI transform was specified, but %s does not look "
                         "pre-transformed. Assuming non-umi data." % fq1)
            return data

    if file_exists(transform):
        transform_file = transform
    else:
        transform_file = get_transform_file(transform)
        if not file_exists(transform_file):
            logger.error(
                "The UMI transform can be specified as either a file or a "
                "bcbio-supported transform. Either the file %s does not exist "
                "or the transform is not supported by bcbio. Supported "
                "transforms are %s."
                % (dd.get_umi_type(data), ", ".join(SUPPORTED_TRANSFORMS)))
            sys.exit(1)
    out_base = dd.get_sample_name(data) + ".umitransformed.fq.gz"
    out_file = os.path.join(umi_dir, out_base)
    if file_exists(out_file):
        data["files"] = [out_file]
        return data
    umis = config_utils.get_program("umis", data, default="umis")
    cores = dd.get_num_cores(data)
    # skip transformation if the file already looks transformed
    with open_fastq(fq1) as in_handle:
        read = next(in_handle)
        if "UMI_" in read:
            data["files"] = [out_file]
            return data

    cmd = ("{umis} fastqtransform {transform_file} "
           "--cores {cores} "
           "{fq1}"
           "| seqtk seq -L 20 - | gzip > {tx_out_file}")
    message = ("Inserting UMI and barcode information into the read name of %s"
               % fq1)
    with file_transaction(out_file) as tx_out_file:
        do.run(cmd.format(**locals()), message)
    data["files"] = [out_file]
    return data
