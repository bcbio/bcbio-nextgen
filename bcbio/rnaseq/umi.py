"""
Unique Molecular Identifier (UMI) handling
Most of this either uses Valentine Svennson's umis repository or adapts
code written from it.
https://github.com/vals/umis
"""
import pandas as pd
import json
import os
import copy
import glob
from itertools import repeat, chain

import bcbio.pipeline.datadict as dd
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.utils import (file_exists, safe_makedir, is_gzipped)
from bcbio.distributed.transaction import file_transaction
from bcbio.bam.fastq import open_fastq
from bcbio.log import logger

transforms = {"harvard-indrop":
              {"read1": r"""(?P<name>^@.*)\n(?P<CB1>\w{8,11})(GAGTGATTGCTTGTGACGCCTT){s<=3}(?P<CB2>\w{8})(?P<MB>\w{6})(.*)\n+(.*)\n(.*)\n""",
               "read2": r"""(@.*)\n(?P<seq>.*)\n\+(.*)\n(?P<qual>.*)\n"""},
              "harvard-indrop-v2":
              {"read2": r"""(?P<name>^@.*)\n(?P<CB1>\w{8,11})(GAGTGATTGCTTGTGACGCCTT){s<=3}(?P<CB2>\w{8})(?P<MB>\w{6})(.*)\n+(.*)\n(.*)\n""",
                  "read1": r"""(@.*)\n(?P<seq>.*)\n\+(.*)\n(?P<qual>.*)\n"""},
              "CEL-seq":
              {"read1": r"""(?P<name>@.*) .*\n(?P<CB>.{8})(?P<MB>.{4})(.*)\n\+(.*)\n(.*)\n""",
               "read2": r"""(@.*)\n(?P<seq>.*)\n\+(.*)\n(?P<qual>.*)\n"""},
              "harvard-indrop-v3":
              {"read1": r"""(?P<name>[^\s]+).*\n(?P<seq>.*)\n\+(.*)\n(?P<qual>.*)\n""",
               "read2": r"""(.*)\n(?P<CB1>.*)\n(.*)\n(.*)\n""",
               "read3": r"""(.*)\n(?P<SB>.*)\n(.*)\n(.*)\n""",
               "read4": r"""(.*)\n(?P<CB2>.{8})(?P<MB>.{6})\n(.*)\n(.*)\n"""}}

def write_transform_file(transform_data, out_file):
    """
    write out the regex to pull out the UMI and cellular barcodes from
    the reads to a JSON file, for use with umis.py
    """
    if file_exists(out_file):
        return out_file

    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            json.dump(transform_data, out_handle)
    return out_file

def umi_transform(data):
    """
    transform each read by identifying the barcode and UMI for each read
    and putting the information in the read name
    """
    fqfiles = data["files"]
    fqfiles.extend(list(repeat("", 4-len(fqfiles))))
    fq1, fq2, fq3, fq4 = fqfiles
    umi_dir = os.path.join(dd.get_work_dir(data), "umis")
    safe_makedir(umi_dir)
    transform = dd.get_umi_type(data)
    if file_exists(transform):
        transform_file = transform
    else:
        transform_data = transforms[transform]
        transform_file = os.path.join(umi_dir, transform + ".json")
        transform_file = write_transform_file(transform_data, transform_file)
    out_base = dd.get_sample_name(data) + ".umitransformed.fq.gz"
    out_file = os.path.join(umi_dir, out_base)
    if file_exists(out_file):
        data["files"] = [out_file]
        return [[data]]
    if len(dd.get_cellular_barcodes(data)) == 2:
        split_option = "--separate_cb"
    else:
        split_option = ""
    umis = config_utils.get_program("umis", data, default="umis")
    cores = dd.get_num_cores(data)
    # skip transformation if the file already looks transformed
    with open_fastq(fq1) as in_handle:
        read = in_handle.next()
        if "UMI_" in read:
            data["files"] = [out_file]
            return [[data]]

    cmd = ("{umis} fastqtransform {split_option} {transform_file} "
           "--cores {cores} "
           "{fq1} {fq2} {fq3} {fq4}"
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
    correction = dd.get_cellular_barcode_correction(data)
    bc = dd.get_cellular_barcodes(data)
    if not bc:
        return [[data]]
    bc1 = None
    bc2 = None
    umi_dir = os.path.join(dd.get_work_dir(data), "umis")
    if isinstance(bc, basestring):
        bc1 = bc
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
        if correction:
            cmd += "--nedit {correction} "
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
    fq1_cmd = fq1
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
    positional = "--positional" if dd.get_positional_umi(data, False) else ""
    message = "Counting alignments of transcripts in %s." % bam
    cmd = ("{umis} tagcount {positional} --cb_cutoff {cutoff} --cb_histogram "
           "{cb_histogram} {bam} {tx_out_file}")
    with file_transaction(out_file) as tx_out_file:
        do.run(cmd.format(**locals()), message)
    data = dd.set_count_file(data, out_file)
    return [[data]]

def get_barcode_metadata(data):
    barcode_file = dd.get_barcode_file(data)
    df = pd.read_csv(barcode_file, sep=",", header=0)
    barcodes = df["barcodes"]

def convert_to_kallisto(data):
    files = dd.get_input_sequence_files(data)
    if len(files) == 2:
        fq1, fq2 = files
    else:
        fq1, fq2 = files[0], None
    samplename = dd.get_sample_name(data)
    work_dir = dd.get_work_dir(data)
    kallisto_dir = os.path.join(work_dir, "kallisto", samplename, "fastq")
    out_file = os.path.join(kallisto_dir, "barcodes.batch")
    umis = config_utils.get_program("umis", dd.get_config(data))
    if file_exists(out_file):
        return out_file
    if dd.get_minimum_barcode_depth(data):
        cb_histogram = os.path.join(work_dir, "umis", samplename, "cb-histogram.txt")
        cb_cutoff = dd.get_minimum_barcode_depth(data)
        cb_options = "--cb_histogram {cb_histogram} --cb_cutoff {cb_cutoff}"
        cb_options = cb_options.format(**locals())
    else:
        cb_options = ""
    cmd = ("{umis} kallisto {cb_options} --out_dir {tx_kallisto_dir} {fq1}")
    with file_transaction(data, kallisto_dir) as tx_kallisto_dir:
        safe_makedir(tx_kallisto_dir)
        message = ("Transforming %s to Kallisto singlecell format. "
                   % fq1)
        do.run(cmd.format(**locals()), message)
    return out_file

def demultiplex_samples(data):
    """
    demultiplex a fastqtransformed FASTQ file into separate sample barcode files
    """
    files = data["files"]
    if len(files) == 2:
        logger.error("Sample demultiplexing doesn't handle paired-end reads, but "
            "we can add it. Open an issue here https://github.com/chapmanb/bcbio-nextgen/issues if you need this and we'll add it.")
        sys.exit(1)
    else:
        fq1 = files[0]
    # check if samples need to be demultiplexed
    with open_fastq(fq1) as in_handle:
        read = in_handle.next()
        if "SAMPLE_" not in read:
            return [[data]]
    bcfile = dd.get_sample_barcodes(data)
    if not bcfile:
        logger.error("Sample demultiplexing needs a list of known indexes provided "
                     "with via the sample_barcodes option in the algorithm section.")
        sys.exit(1)
    work_dir = os.path.join(dd.get_work_dir(data), "umis")
    sample_dir = os.path.join(work_dir, dd.get_sample_name(data))
    demulti_dir = os.path.join(sample_dir, "demultiplexed")
    demultiplexed = glob.glob(os.path.join(demulti_dir, "*.fq*"))
    if demultiplexed:
        return [split_demultiplexed_sampledata(data, demultiplexed)]
    umis = config_utils.get_program("umis", data, default="umis")
    cmd = ("{umis} demultiplex_samples --nedit 1 --barcodes {bcfile} "
           "--out_dir {tx_dir} {fq1}")
    msg = "Demultiplexing {fq1}."
    with file_transaction(data, demulti_dir) as tx_dir:
        do.run(cmd.format(**locals()), msg.format(**locals()))
    demultiplexed = glob.glob(os.path.join(demulti_dir, "*.fq*"))
    return [split_demultiplexed_sampledata(data, demultiplexed)]

def split_demultiplexed_sampledata(data, demultiplexed):
    """
    splits demultiplexed samples into separate entries in the global sample
    datadict
    """
    datadicts = []
    samplename = dd.get_sample_name(data)
    for fastq in demultiplexed:
        barcode = os.path.basename(fastq).split(".")[0]
        datadict = copy.deepcopy(data)
        datadict = dd.set_sample_name(datadict, samplename + "-" + barcode)
        datadict = dd.set_description(datadict, samplename + "-" + barcode)
        datadict["rgnames"]["rg"] = samplename + "-" + barcode
        datadict["name"]= ["", samplename + "-" + barcode]
        datadict["files"] = [fastq]
        datadicts.append(datadict)
    return datadicts
