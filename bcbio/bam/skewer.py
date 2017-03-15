import os
import shutil
import bcbio.pipeline.datadict as dd
from bcbio.pipeline import config_utils
from bcbio.distributed.transaction import tx_tmpdir
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
    mode = "tail" if not fq2 else "pe"
    cmd = ("{skewer} --min 25 --threads {nthreads} "
           "{fw_cmd} "
           "{rv_cmd} "
           "-m {mode} "
           "--compress --output {out_stem} ") + file_string
    with tx_tmpdir(data, out_dir) as tx_out_dir:
        safe_makedir(tx_out_dir)
        out_stem = os.path.join(tx_out_dir, samplename)
        fw_cmd = _fw_command(data, tx_out_dir)
        rv_cmd = _rv_command(data, tx_out_dir)
        message = "Trimming {fq1}, {fq2} with skewer.".format(**locals())
        do.run(cmd.format(**locals()), message)
        shutil.move(os.path.join(tx_out_dir, os.path.basename(of1)), out_dir)
        shutil.move(os.path.join(tx_out_dir, os.path.basename(of2)), out_dir)
        shutil.move(out_stem + "-trimmed.log", out_dir)
    return of1, of2

def _fw_command(data, tx_out_dir):
    fw = _get_fw_adapters(data)
    if not fw:
        return ""
    else:
        fw_file = _write_fasta_file(fw, "forward", tx_out_dir)
        return "-x %s" % fw_file

def _rv_command(data, tx_out_dir):
    rv = _get_rv_adapters(data)
    if not rv:
        return ""
    else:
        rv_file = _write_fasta_file(rv, "reverse", tx_out_dir)
        return "-y %s" % rv_file

def _write_fasta_file(adapters, ext, tx_out_dir):
    tfile = os.path.join(tx_out_dir, "adapters-%s.fa" % ext)
    with open(tfile, "w") as out_handle:
        for i, adapter in enumerate(adapters):
            if file_exists(adapter):
                with open(adapter) as in_handle:
                    for line in in_handle:
                        out_handle.write(line)
            else:
                out_handle.write(">%d\n" % i)
                out_handle.write(adapter + "\n")
    return tfile

def _get_fw_adapters(data):
    builtin = [FW_ADAPTERS[x] for x in dd.get_adapters(data) if x in FW_ADAPTERS]
    return flatten(builtin + dd.get_custom_trim(data))

def _get_rv_adapters(data):
    builtin = [RV_ADAPTERS[x] for x in dd.get_adapters(data) if x in FW_ADAPTERS]
    return flatten(builtin + dd.get_custom_trim(data))
