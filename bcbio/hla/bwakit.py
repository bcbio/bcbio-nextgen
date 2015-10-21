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
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do

def run(data):
    """HLA typing with bwakit, parsing output from called genotype files.
    """
    bwakit_dir = os.path.dirname(os.path.realpath(utils.which("run-bwamem")))
    align_file = dd.get_align_bam(data)
    hla_base = os.path.join(utils.safe_makedir(os.path.join(os.path.dirname(align_file), "hla")),
                            os.path.basename(align_file) + ".hla")
    if len(glob.glob(hla_base + ".*")) > 0:
        out_file = hla_base + ".top"
        if not utils.file_exists(out_file):
            cmd = "{bwakit_dir}/run-HLA {hla_base}"
            #do.run(cmd.format(**locals()), "HLA typing with bwakit")
            out_file = _organize_calls(out_file, hla_base, data)
        data["hla"] = {"calls": out_file}
    return data

def _organize_calls(out_file, hla_base, data):
    """Prepare genotype calls, reporting best call along with quality metrics.
    """
    hla_truth = get_hla_truthset(data)
    align_file = dd.get_align_bam(data)
    sample = dd.get_sample_name(data)
    with file_transaction(data, out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            writer = csv.writer(out_handle)
            writer.writerow(["sample", "locus", "mismatches", "options", "alleles", "expected"])
            for genotype_file in glob.glob("%s.HLA-*.gt" % (hla_base)):
                hla_locus = os.path.basename(genotype_file).replace(
                        "%s.hla.HLA-" % os.path.basename(align_file), "").replace(".gt", "")
                with open(genotype_file) as in_handle:
                    total_options = 0
                    for i, line in enumerate(in_handle):
                        if i == 0:
                            _, allele_one, allele_two, mismatches = line.split("\t")[:4]
                        total_options += 1
                    if total_options > 0:
                        truth_alleles = tz.get_in([sample, hla_locus], hla_truth, [])
                        writer.writerow([sample, hla_locus, mismatches, total_options,
                                         ";".join([allele_one, allele_two]),
                                         ";".join(truth_alleles)])
    return out_file

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
                out = tz.update_in(out, [sample, locus], lambda x: alleles.split(";"))
    return out
