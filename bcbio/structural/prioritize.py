"""Prioritize structural variants based on biological information.

Provides high level summaries of structural variants in regions of interest,
as defined by the input configuration. Tries to narrow structural variant calls
based on potential biological targets.
"""
import os
import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import bedutils, vcfutils
from bcbio.structural import lumpy, regions

POST_PRIOR_FNS = {"lumpy": lumpy.run_svtyper_prioritize}

def run(items):
    assert len(items) == 1, "Expect one input to biological prioritization"
    data = items[0]
    inputs = []
    for call in data.get("sv", []):
        vcf_file = call.get("vcf_file", call.get("vrn_file"))
        if vcf_file and vcf_file.endswith((".vcf", "vcf.gz")):
            pp_fn = POST_PRIOR_FNS.get(call["variantcaller"])
            if pp_fn:
                pp_fn = pp_fn(call)
            inputs.append((call["variantcaller"], vcf_file, pp_fn))
    if len(inputs) > 0:
        prioritize_by = tz.get_in(["config", "algorithm", "svprioritize"], data)
        if not prioritize_by:
            raise ValueError("Missing structural variant prioritization with `svprioritize`")
        work_dir = _sv_workdir(data)
        priority_files = [_prioritize_vcf(vcaller, vfile, prioritize_by, post_prior_fn, work_dir, data)
                          for vcaller, vfile, post_prior_fn in inputs]
        priority_tsv = _combine_files(priority_files, work_dir, data)
        data["sv"].append({"variantcaller": "sv-prioritize", "vrn_file": priority_tsv})
    return [data]

def _prioritize_vcf(caller, vcf_file, prioritize_by, post_prior_fn, work_dir, data):
    """Provide prioritized tab delimited output for a single caller.
    """
    sample = dd.get_sample_name(data)
    out_file = os.path.join(work_dir, "%s-%s-prioritize.tsv" % (sample, caller))
    if not utils.file_exists(out_file):
        priority_vcf = "%s.vcf.gz" % utils.splitext_plus(out_file)[0]
        if not utils.file_exists(priority_vcf):
            with file_transaction(data, priority_vcf) as tx_out_file:
                resources = config_utils.get_resources("bcbio_prioritize", data["config"])
                jvm_opts = " ".join(resources.get("jvm_opts", ["-Xms1g", "-Xmx4g"]))
                cmd = ("bcbio-prioritize {jvm_opts} known -i {vcf_file} -o {tx_out_file} -k {prioritize_by}")
                do.run(cmd.format(**locals()), "Prioritize: select in known regions of interest")
        if post_prior_fn:
            priority_vcf = post_prior_fn(priority_vcf, work_dir, data)
        simple_vcf = "%s-simple.vcf.gz" % utils.splitext_plus(priority_vcf)[0]
        if not utils.file_exists(simple_vcf):
            with file_transaction(data, simple_vcf) as tx_out_file:
                transcript_file = regions.get_sv_bed(data, "transcripts1000", work_dir)
                if transcript_file:
                    transcript_file = vcfutils.bgzip_and_index(transcript_file, data["config"])
                    ann_opt = "--gene_bed %s" % transcript_file
                else:
                    ann_opt = ""
                cmd = "simple_sv_annotation.py {ann_opt} -o - {priority_vcf} | bgzip -c > {tx_out_file}"
                do.run(cmd.format(**locals()), "Prioritize: simplified annotation output")
        simple_vcf = vcfutils.bgzip_and_index(vcfutils.sort_by_ref(simple_vcf, data), data["config"])
        with file_transaction(data, out_file) as tx_out_file:
            cmd = ("zcat {simple_vcf} | vawk -v SNAME={sample} -v CALLER={caller} "
                   """'{{if (($7 == "PASS" || $7 == ".") && (S${sample}$GT != "0/0")) """
                   "print CALLER,SNAME,$1,$2,I$END,"
                   """I$SVTYPE=="BND" ? I$SVTYPE":"$3":"I$MATEID : I$SVTYPE,"""
                   "I$KNOWN,I$END_GENE,I$LOF,I$SIMPLE_ANN,"
                   "S${sample}$SR,S${sample}$PE}}' > {tx_out_file}")
            do.run(cmd.format(**locals()), "Prioritize: convert to tab delimited")
    return out_file

def _combine_files(tsv_files, work_dir, data):
    """Combine multiple priority tsv files into a final sorted output.
    """
    header = "\t".join(["caller", "sample", "chrom", "start", "end", "svtype", "known", "end_gene",
                        "lof", "annotation", "split_read_support", "paired_end_support"])
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
