import os
import tempfile
import bcbio.pipeline.datadict as dd
from bcbio.pipeline import config_utils
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.utils import safe_makedir, file_exists, flatten

FW_ADAPTERS = {"truseq": ["AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"]}
RV_ADAPTERS = {"truseq": ["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"]}

def trim_adapters(data):
    fq1, fq2 = dd.get_input_sequence_files(data)
    skewer = config_utils.get_program("skewer", data, default="skewer")
    nthreads = dd.get_num_cores(data)
    samplename = dd.get_sample_name(data)
    out_dir = safe_makedir(os.path.join(dd.get_work_dir(data), "trimmed", samplename))
    of1 = os.path.join(out_dir, samplename + "-trimmed-pair1.fastq.gz")
    of2 = os.path.join(out_dir, samplename + "-trimmed-pair2.fastq.gz")
    of2 = of2 if fq2 else None
    if fq1 and fq2:
        if file_exists(of1) and file_exists(of2):
            return of1, of2
    else:
        if file_exists(of1):
            return of1, None
    safe_makedir(out_dir)
    file_string = "{fq1} {fq2} " if fq2 else "{fq1} "
    fw_cmd = _fw_command(data)
    rv_cmd = _rv_command(data)
    mode = "tail" if not fq2 else "pe"
    cmd = ("{skewer} --min 25 --threads {nthreads} -q 5 "
           "{fw_cmd} "
           "{rv_cmd} "
           "-m {mode} "
           "--compress --output {out_stem} ") + file_string
    with file_transaction(out_dir) as tx_out_dir:
        safe_makedir(tx_out_dir)
        out_stem = os.path.join(tx_out_dir, samplename)
        message = "Trimming {fq1}, {fq2} with skewer.".format(**locals())
        do.run(cmd.format(**locals()), message)
    return of1, of2

def _fw_command(data):
    fw = _get_fw_adapters(data)
    if not fw:
        return ""
    else:
        fw_file = _write_fasta_file(fw)
        return "-x %s" % fw_file

def _rv_command(data):
    rv = _get_rv_adapters(data)
    if not rv:
        return ""
    else:
        rv_file = _write_fasta_file(rv)
        return "-y %s" % rv_file

def _write_fasta_file(adapters):
    count = 0
    tfile = tempfile.NamedTemporaryFile(delete=False).name
    with open(tfile, "w") as out_handle:
        for adapter in adapters:
            out_handle.write(">%d\n" % count)
            out_handle.write(adapter + "\n")
            count += 1
    return tfile

def _get_fw_adapters(data):
    builtin = [FW_ADAPTERS[x] for x in dd.get_adapters(data) if x in FW_ADAPTERS]
    return flatten(builtin + dd.get_custom_trim(data))

def _get_rv_adapters(data):
    builtin = [RV_ADAPTERS[x] for x in dd.get_adapters(data) if x in FW_ADAPTERS]
    return flatten(builtin + dd.get_custom_trim(data))
