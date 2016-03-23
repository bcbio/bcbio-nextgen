"""
Unique Molecular Identifier (UMI) handling
Most of this either uses Valentine Svennson's umis repository or adapts
code written from it.
https://github.com/vals/umis
"""
import json
import os

import bcbio.pipeline.datadict as dd
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.utils import (file_exists, safe_makedir, is_gzipped)
from bcbio.distributed.transaction import file_transaction


transforms = {"harvard-indrop":
              {"json": {
                  "read1": r"""(?P<name>^@.*)\n(?P<CB1>\w{8,11})(GAGTGATTGCTTGTGACGCCTT){s<=3}(?P<CB2>\w{8})(?P<MB>\w{6})(.*)\n+(.*)\n(.*)\n""",
                  "read2": r"""(@.*)\n(?P<seq>.*)\n\+(.*)\n(?P<qual>.*)\n"""},
               "dual": True},
              "harvard-indrop-v2":
              {"json": {
                  "read2": r"""(?P<name>^@.*)\n(?P<CB1>\w{8,11})(GAGTGATTGCTTGTGACGCCTT){s<=3}(?P<CB2>\w{8})(?P<MB>\w{6})(.*)\n+(.*)\n(.*)\n""",
                  "read1": r"""(@.*)\n(?P<seq>.*)\n\+(.*)\n(?P<qual>.*)\n"""},
               "dual": True},
              "CEL-seq":
              {"json": {
                  "read1": r"""(?P<name>@.*) .*\\n(?P<CB>.{8})(?P<MB>.{4})(.*)\\n\\+(.*)\\n(.*)\\n""",
                  "read2": r"""(@.*)\\n(?P<seq>.*)\\n\\+(.*)\\n(?P<qual>.*)\\n"""},
               "dual": False}}

def write_transform_file(transform_data, out_file):
    """
    write out the regex to pull out the UMI and cellular barcodes from
    the reads to a JSON file, for use with umis.py
    """
    if file_exists(out_file):
        return out_file

    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            json.dump(transform_data["json"], out_handle)
    return out_file

def umi_transform(data):
    """
    transform each read by identifying the barcode and UMI for each read
    and putting the information in the read name
    """
    fq1, fq2 = dd.get_input_sequence_files(data)
    fq2 = fq2 if fq2 else ""
    umi_dir = os.path.join(dd.get_work_dir(data), "umis")
    transform = dd.get_umi_type(data)
    transform_data = transforms[transform]
    safe_makedir(umi_dir)
    transform_file = os.path.join(umi_dir, transform + ".json")
    transform_file = write_transform_file(transform_data, transform_file)
    out_base = dd.get_sample_name(data) + ".umitransformed.fq.gz"
    out_file = os.path.join(umi_dir, out_base)
    if file_exists(out_file):
        data["files"] = [out_file]
        return [[data]]
    index_option = "--dual_index" if transform_data["dual"] else ""
    if len(dd.get_cellular_barcodes(data)) == 2:
        split_option = "--separate_cb"
    else:
        split_option = ""
    umis = config_utils.get_program("umis", data, default="umis")
    cmd = ("{umis} fastqtransform {index_option} {split_option} {transform_file} "
           "{fq1} {fq2} "
           "| seqtk seq -L 20 - | gzip > {tx_out_file}")
    message = ("Inserting UMI and barcode information into the read name of %s"
               % fq1)
    with file_transaction(out_file) as tx_out_file:
        do.run(cmd.format(**locals()), message)
    data["files"] = [out_file]
    return [[data]]

def filter_barcodes(data):
    fq1 = dd.get_input_sequence_files(data)[0]
    umi_dir = os.path.join(dd.get_work_dir(data), "umis")
    bc = dd.get_cellular_barcodes(data)
    if not bc:
        return [[data]]
    bc1 = None
    bc2 = None
    umi_dir = os.path.join(dd.get_work_dir(data), "umis")
    if len(bc) == 1:
        bc1 = bc[0]
    if len(bc) == 2:
        bc1 = bc[0]
        bc2 = bc[1]
    out_base = dd.get_sample_name(data) + ".filtered.fq.gz"
    out_file = os.path.join(umi_dir, out_base)
    if file_exists(out_file):
        data["files"] = [out_file]
        return [[data]]

    ncores = dd.get_num_cores(data)
    cmd = "{umis} cb_filter --cores {ncores} "
    if bc1:
        cmd += "--bc1 {bc1} "
    if bc2:
        cmd += "--bc2 {bc2} "

    fq1_cmd = "{fq1} " if not is_gzipped(fq1) else "<(gzip -cd {fq1}) "
    fq1_cmd = fq1_cmd.format(fq1=fq1)
    cmd += "{fq1_cmd} | gzip > {tx_out_file}"

    sample_dir = os.path.join(umi_dir, dd.get_sample_name(data))
    safe_makedir(sample_dir)
    umis = config_utils.get_program("umis", data, default="umis")
    with file_transaction(out_file) as tx_out_file:
        message = "Filtering by cellular barcode."
        do.run(cmd.format(**locals()), message)
    data["files"] = [out_file]
    return [[data]]

def barcode_histogram(data):
    fq1 = dd.get_input_sequence_files(data)[0]
    umi_dir = os.path.join(dd.get_work_dir(data), "umis")
    sample_dir = os.path.join(umi_dir, dd.get_sample_name(data))
    umis = config_utils.get_program("umis", data, default="umis")
    safe_makedir(sample_dir)
    out_file = os.path.join(sample_dir, "cb-histogram.txt")
    if file_exists(out_file):
        return [[data]]
    fq1_cmd = "{fq1}" if not is_gzipped(fq1) else "<(gzip -cd {fq1})"
    fq1_cmd = fq1_cmd.format(fq1=fq1)
    cmd = "{umis} cb_histogram {fq1_cmd} > {tx_out_file}"
    with file_transaction(out_file) as tx_out_file:
        message = "Computing cellular barcode counts for %s." % fq1
        do.run(cmd.format(**locals()), message)
    return [[data]]

def tagcount(data):
    bam = dd.get_transcriptome_bam(data)
    umi_dir = os.path.join(dd.get_work_dir(data), "umis")
    sample_dir = os.path.join(umi_dir, dd.get_sample_name(data))
    out_file = os.path.join(sample_dir, dd.get_sample_name(data) + ".counts")
    if file_exists(out_file):
        data = dd.set_count_file(data, out_file)
        return [[data]]
    umis = config_utils.get_program("umis", data, default="umis")
    safe_makedir(sample_dir)
    cutoff = dd.get_minimum_barcode_depth(data)
    cb_histogram = os.path.join(sample_dir, "cb-histogram.txt")
    message = "Counting alignments of transcripts in %s." % bam
    cmd = ("{umis} tagcount --positional --cb_cutoff {cutoff} --cb_histogram "
           "{cb_histogram} {bam} {tx_out_file}")
    with file_transaction(out_file) as tx_out_file:
        do.run(cmd.format(**locals()), message)
    data = dd.set_count_file(data, out_file)
    return [[data]]

def get_barcode_metadata(data):
    barcode_file = dd.get_barcode_file(data)
    df = pd.read_csv(barcode_file, sep=",", header=0)
    barcodes = df["barcodes"]
