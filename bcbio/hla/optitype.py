"""HLA typing from extracted HLA alleles with OptiType

https://github.com/FRED-2/OptiType
"""
import csv
import glob
import os
import re
import sys
import shutil

import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.hla import bwakit
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do

SUPPORTED_HLAS = ["HLA-A", "HLA-B", "HLA-C"]

def run(data):
    """HLA typing with OptiType, parsing output from called genotype files.
    """
    hlas = []
    for hla_fq in tz.get_in(["hla", "fastq"], data, []):
        hla_type = re.search(r"[.-](?P<hlatype>HLA-[\w-]+).fq", hla_fq).group("hlatype")
        if hla_type in SUPPORTED_HLAS:
            if utils.file_exists(hla_fq):
                hlas.append((hla_type, hla_fq))
    if len(hlas) > 0:
        out_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "align",
                                                  dd.get_sample_name(data), "hla",
                                                  "OptiType-HLA-A_B_C"))
        # When running UMIs and hla typing we want to pick the original fastqs
        if len(hlas) > len(SUPPORTED_HLAS):
            hlas = [x for x in hlas if os.path.basename(x[1]).find("-cumi") == -1]
        if len(hlas) == len(SUPPORTED_HLAS):
            hla_fq = combine_hla_fqs(hlas, out_dir + "-input.fq", data)
            if utils.file_exists(hla_fq):
                out_file = glob.glob(os.path.join(out_dir, "*", "*_result.tsv"))
                if len(out_file) > 0:
                    out_file = out_file[0]
                else:
                    out_file = _call_hla(hla_fq, out_dir, data)
                out_file = _prepare_calls(out_file, os.path.dirname(out_dir), data)
                data["hla"].update({"call_file": out_file,
                                    "hlacaller": "optitype"})
    return data

def combine_hla_fqs(hlas, out_file, data):
    """OptiType performs best on a combination of all extracted HLAs.
    """
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                for hla_type, hla_fq in hlas:
                    if utils.file_exists(hla_fq):
                        with open(hla_fq) as in_handle:
                            shutil.copyfileobj(in_handle, out_handle)
    return out_file

def _prepare_calls(result_file, out_dir, data):
    """Write summary file of results of HLA typing by allele.
    """
    sample = dd.get_sample_name(data)
    out_file = os.path.join(out_dir, "%s-optitype.csv" % (sample))
    if not utils.file_uptodate(out_file, result_file):
        hla_truth = bwakit.get_hla_truthset(data)
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                writer = csv.writer(out_handle)
                allele_info = _parse_result_file(result_file)
                if len(allele_info) == 1:
                    writer.writerow(["sample", "locus", "alleles", "expected", "validates"])
                else:
                    writer.writerow(["sample", "local", "index", "alleles", "score"])
                for j, (alleles, score) in enumerate(allele_info):
                    for hla_locus, call_alleles in alleles:
                        truth_alleles = tz.get_in([sample, hla_locus], hla_truth, [])
                        if len(allele_info) == 1:
                            writer.writerow([sample, hla_locus,
                                             ";".join(call_alleles), ";".join(truth_alleles),
                                             bwakit.matches_truth(call_alleles, truth_alleles, data)])
                        else:
                            writer.writerow([sample, hla_locus, j, ";".join(call_alleles), score])
    return out_file

def _parse_result_file(result_file):
    with open(result_file) as in_handle:
        header = in_handle.readline().rstrip("\r\n").split("\t")
        all_locus_is = []
        for hla_locus in SUPPORTED_HLAS:
            hla_locus = hla_locus.replace("HLA-", "")
            all_locus_is.append((hla_locus, [i for i, h in enumerate(header) if h.startswith(hla_locus)]))
        score_i = [i for i, h in enumerate(header) if h == "Objective"][0]
        hits = []
        for line in in_handle:
            hit = line.rstrip("\r\n").split("\t")
            assert len(hit) == len(header), result_file
            cur_hlas = []
            for hla_locus, locus_is in all_locus_is:
                cur_hlas.append((hla_locus, ["HLA-%s" % hit[i] for i in locus_is]))
            hits.append((cur_hlas, hit[score_i]))
        return hits

def _call_hla(hla_fq, out_dir, data):
    """Run OptiType HLA calling for a specific fastq input.
    """
    bin_dir = os.path.dirname(os.path.realpath(sys.executable))
    out_dir = utils.safe_makedir(out_dir)
    with tx_tmpdir(data, os.path.dirname(out_dir)) as tx_out_dir:
        config_file = os.path.join(tx_out_dir, "config.ini")
        with open(config_file, "w") as out_handle:
            razers3 = os.path.join(bin_dir, "razers3")
            if not os.path.exists(razers3):
                raise ValueError("Could not find razers3 executable at %s" % (razers3))
            out_handle.write(CONFIG_TMPL.format(razers3=razers3, cores=dd.get_cores(data)))
        resources = config_utils.get_resources("optitype", data["config"])
        if resources.get("options"):
            opts = " ".join([str(x) for x in resources["options"]])
        else:
            opts = ""
        cmd = ("OptiTypePipeline.py -v --dna {opts} -o {tx_out_dir} "
                "-i {hla_fq} -c {config_file}")
        do.run(cmd.format(**locals()), "HLA typing with OptiType")
        for outf in os.listdir(tx_out_dir):
            shutil.move(os.path.join(tx_out_dir, outf), os.path.join(out_dir, outf))
    out_file = glob.glob(os.path.join(out_dir, "*", "*_result.tsv"))
    assert len(out_file) == 1, "Expected one result file for OptiType, found %s" % out_file
    return out_file[0]


CONFIG_TMPL = """
[mapping]
razers3={razers3}
threads={cores}
[ilp]
solver=cbc
threads={cores}
[behavior]
deletebam=true
unpaired_weight=0
use_discordant=false
"""
