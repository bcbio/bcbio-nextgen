"""Prioritize structural variants based on biological information.

Provides high level summaries of structural variants in regions of interest,
as defined by the input configuration. Tries to narrow structural variant calls
based on potential biological targets.
"""
import os
import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import bedutils, vcfutils

def run(items):
    assert len(items) == 1, "Expect one input to biological prioritization"
    data = items[0]
    inputs = []
    for call in data.get("sv", []):
        vcf_file = call.get("vcf_file", call.get("vrn_file"))
        if vcf_file and vcf_file.endswith((".vcf", "vcf.gz")):
            inputs.append((call["variantcaller"], vcf_file))
    if len(inputs) > 0:
        prioritize_by = tz.get_in(["config", "algorithm", "svprioritize"], data)
        if not prioritize_by:
            raise ValueError("Missing structural variant prioritization with `svprioritize`")
        work_dir = _sv_workdir(data)
        priority_files = [_prioritize_vcf(vcaller, vfile, prioritize_by, work_dir, data)
                          for vcaller, vfile in inputs]
        priority_tsv = _combine_files(priority_files, work_dir, data)
        data["sv"].append({"variantcaller": "sv-prioritize", "vrn_file": priority_tsv})
    return [data]

def _prioritize_vcf(caller, vcf_file, prioritize_by, work_dir, data):
    """Provide prioritized tab delimited output for a single caller.
    """
    sample = dd.get_sample_name(data)
    out_file = os.path.join(work_dir, "%s-%s-prioritize.tsv" % (sample, caller))
    if not utils.file_exists(out_file):
        priority_vcf = "%s.vcf.gz" % utils.splitext_plus(out_file)[0]
        if not utils.file_exists(priority_vcf):
            with file_transaction(data, priority_vcf) as tx_out_file:
                cmd = ("bcbio-prioritize known -i {vcf_file} -o {tx_out_file} -k {prioritize_by}")
                do.run(cmd.format(**locals()), "Prioritize: select in known regions of interest")
        simple_vcf = "%s-simple.vcf.gz" % utils.splitext_plus(priority_vcf)[0]
        if not utils.file_exists(simple_vcf):
            with file_transaction(data, simple_vcf) as tx_out_file:
                cmd = "simple_sv_annotation.py -o - {priority_vcf} | bgzip -c > {tx_out_file}"
                do.run(cmd.format(**locals()), "Prioritize: simplified annotation output")
        simple_vcf = vcfutils.bgzip_and_index(simple_vcf, data["config"])
        with file_transaction(data, out_file) as tx_out_file:
            cmd = ("zcat {simple_vcf} | vawk -v SNAME={sample} -v CALLER={caller} "
                   """'{{if (($7 == "PASS" || $7 == ".")) """
                   "print CALLER,SNAME,$1,$2,I$END,I$SVTYPE,I$KNOWN,I$LOF,I$SIMPLE_ANN,"
                   "S${sample}$SR,S${sample}$PE}}' > {tx_out_file}")
            do.run(cmd.format(**locals()), "Prioritize: convert to tab delimited")
    return out_file

def _combine_files(tsv_files, work_dir, data):
    """Combine multiple priority tsv files into a final sorted output.
    """
    header = "\t".join(["caller", "sample", "chrom", "start", "end", "svtype", "known", "lof",
                        "annotation", "split_read_support", "paired_end_support"])
    sample = dd.get_sample_name(data)
    out_file = os.path.join(work_dir, "%s-prioritize.tsv" % (sample))
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            input_files = " ".join(tsv_files)
            sort_cmd = bedutils.get_sort_cmd()
            cmd = "{{ echo '{header}'; cat {input_files} | {sort_cmd} -k3,3 -k4,4n; }} > {tx_out_file}"
            do.run(cmd.format(**locals()), "Combine prioritized from multiple callers")
    return out_file

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "prioritize"))
