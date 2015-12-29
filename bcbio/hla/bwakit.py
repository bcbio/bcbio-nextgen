"""Call HLA alleles with assembly methods implemented in bwakit.

https://github.com/lh3/bwa/blob/master/README-alt.md#hla-typing
https://github.com/lh3/bwa/tree/master/bwakit
"""
import csv
import glob
import os

import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.hla import groups as hla_groups
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do

def run(data):
    """HLA typing with bwakit, parsing output from called genotype files.
    """
    bwakit_dir = os.path.dirname(os.path.realpath(utils.which("run-bwamem")))
    hla_fqs = tz.get_in(["hla", "fastq"], data, [])
    if len(hla_fqs) > 0:
        hla_base = os.path.commonprefix(hla_fqs)
        while hla_base.endswith("."):
            hla_base = hla_base[:-1]
        out_file = hla_base + ".top"
        if not utils.file_exists(out_file):
            cmd = "{bwakit_dir}/run-HLA {hla_base}"
            do.run(cmd.format(**locals()), "HLA typing with bwakit")
            out_file = _organize_calls(out_file, hla_base, data)
        data["hla"].update({"call_file": out_file,
                            "hlacaller": "bwakit"})
    return data

def _organize_calls(out_file, hla_base, data):
    """Prepare genotype calls, reporting best call along with quality metrics.
    """
    hla_truth = get_hla_truthset(data)
    sample = dd.get_sample_name(data)
    with file_transaction(data, out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            writer = csv.writer(out_handle)
            writer.writerow(["sample", "locus", "mismatches", "options", "alleles", "p-groups", "expected",
                             "validates"])
            for genotype_file in glob.glob("%s.HLA-*.gt" % (hla_base)):
                hla_locus = os.path.basename(genotype_file).replace(
                        "%s.HLA-" % os.path.basename(hla_base), "").replace(".gt", "")
                with open(genotype_file) as in_handle:
                    total_options = set([])
                    for i, line in enumerate(in_handle):
                        _, aone, atwo, m = line.split("\t")[:4]
                        pgroups = (hla_groups.hla_protein(aone, data), hla_groups.hla_protein(atwo, data))
                        if i == 0:
                            call_alleles = [aone, atwo]
                            call_pgroups = pgroups
                            mismatches = m
                        total_options.add(pgroups)
                    if len(total_options) > 0:
                        truth_alleles = tz.get_in([sample, hla_locus], hla_truth, [])
                        writer.writerow([sample, hla_locus, mismatches, len(total_options),
                                         ";".join(call_alleles), ";".join(call_pgroups),
                                         ";".join(truth_alleles), matches_truth(call_alleles, truth_alleles, data)])
    return out_file

def matches_truth(call_alleles, truth_alleles, data):
    """Flexibly check if truth and call alleles match, using p-groups.
    """
    if not truth_alleles:
        return ""
    else:
        def _remove_p(x):
            return x[:-1] if x.endswith("P") else x
        t_cmp = set([_remove_p(hla_groups.hla_protein(x, data)) for x in truth_alleles])
        c_cmp = set([_remove_p(hla_groups.hla_protein(x, data)) for x in call_alleles])
        return "yes" if len(t_cmp.intersection(c_cmp)) == len(t_cmp) else "no"

def get_hla_truthset(data):
    """Retrieve expected truth calls for annotating HLA called output.
    """
    val_csv = tz.get_in(["config", "algorithm", "hlavalidate"], data)
    out = {}
    if val_csv and utils.file_exists(val_csv):
        with open(val_csv) as in_handle:
            reader = csv.reader(in_handle)
            reader.next() # header
            for sample, locus, alleles in (l for l in reader if l):
                out = tz.update_in(out, [sample, locus], lambda x: [x.strip() for x in alleles.split(";")])
    return out
