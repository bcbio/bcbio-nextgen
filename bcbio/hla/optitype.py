"""HLA typing from extracted HLA alleles with OptiType

https://github.com/FRED-2/OptiType
"""
import csv
import glob
import os
import sys
import shutil

import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.hla import bwakit
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do

SUPPORTED_HLAS = ["HLA-A", "HLA-B", "HLA-C"]

def run(data):
    """HLA typing with OptiType, parsing output from called genotype files.
    """
    align_file = dd.get_align_bam(data)
    hla_dir = os.path.join(os.path.dirname(align_file), "hla")
    hla_base = os.path.join(hla_dir, os.path.basename(align_file) + ".hla")
    hlas = []
    for hla_fq in glob.glob(hla_base + ".*.fq"):
        hla_type = os.path.splitext(os.path.splitext(os.path.basename(hla_fq))[0])[1].replace(".", "")
        if hla_type in SUPPORTED_HLAS:
            hlas.append((hla_type, hla_fq))
    if len(hlas) > 0:
        hla_calls = []
        for hla_type, hla_fq in hlas:
            out_dir = os.path.join(hla_dir, "OptiType-%s" % hla_type)
            out_file = glob.glob(os.path.join(out_dir, "*", "*_result.tsv"))
            if len(out_file) > 0:
                out_file = out_file[0]
            else:
                out_file = _call_hla(hla_fq, out_dir, data)
            hla_calls.append((hla_type.replace("HLA-", ""), out_file))
        out_file = _combine_calls(hla_calls, hla_dir, data)
        data["hla"] = {"call_file": out_file,
                       "hlacaller": "optitype"}
    return data

def _combine_calls(hla_calls, hla_dir, data):
    """Write summary file of results of HLA typing by allele.
    """
    sample = dd.get_sample_name(data)
    out_file = os.path.join(hla_dir, "%s-optitype.csv" % (sample))
    if not utils.file_uptodate(out_file, hla_calls[0][1]):
        hla_truth = bwakit.get_hla_truthset(data)
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                writer = csv.writer(out_handle)
                writer.writerow(["sample", "locus", "alleles", "expected", "validates"])
                for hla_locus, result_file in hla_calls:
                    truth_alleles = tz.get_in([sample, hla_locus], hla_truth, [])
                    call_alleles, score = _parse_result_file(result_file, hla_locus)
                    writer.writerow([sample, hla_locus,
                                     ";".join(call_alleles), ";".join(truth_alleles),
                                     bwakit.matches_truth(call_alleles, truth_alleles, data)])
    return out_file

def _parse_result_file(result_file, hla_locus):
    with open(result_file) as in_handle:
        header = in_handle.readline().rstrip("\r\n").split("\t")
        locus_is = [i for i, h in enumerate(header) if h.startswith(hla_locus)]
        score_i = [i for i, h in enumerate(header) if h == "Objective"][0]
        hit = in_handle.readline().rstrip("\r\n").split("\t")
        assert len(hit) == len(header), result_file
        return ["HLA-%s" % hit[i] for i in locus_is], hit[score_i]

def _call_hla(hla_fq, out_dir, data):
    """Run OptiType HLA calling for a specific
    """
    bin_dir = os.path.dirname(os.path.realpath(sys.executable))
    with tx_tmpdir(data, os.path.dirname(out_dir)) as tx_out_dir:
        config_file = os.path.join(tx_out_dir, "config.ini")
        with open(config_file, "w") as out_handle:
            razers3 = os.path.join(bin_dir, "razers3")
            if not os.path.exists(razers3):
                raise ValueError("Could not find razers3 executable at %s" % (razers3))
            out_handle.write(CONFIG_TMPL.format(razers3=razers3, cores=dd.get_cores(data)))
        cmd = ("OptiTypePipeline.py -v --dna -o {tx_out_dir} "
                "-i {hla_fq} -c {config_file}")
        do.run(cmd.format(**locals()), "HLA typing with OptiType")
        shutil.move(tx_out_dir, out_dir)
    out_file = glob.glob(os.path.join(out_dir, "*", "*_result.tsv"))
    assert len(out_file) == 1, "Expected one result file for OptiType, found %s" % out_file
    return out_file[0]

CONFIG_TMPL = """
[mapping]
razers3={razers3}
threads={cores}
[ilp]
solver=glpk
threads=1
[behavior]
deletebam=true
unpaired_weight=0
use_discordant=false
"""
