"""Prioritize structural variants based on biological information.

Provides high level summaries of structural variants in regions of interest,
as defined by the input configuration. Tries to narrow structural variant calls
based on potential biological targets.
"""
import os
import pandas as pd
import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import bedutils, vcfutils

POST_PRIOR_FNS = {}

def run(items):
    assert len(items) == 1, ("Expect one input to biological prioritization: %s" %
                             ", ".join([dd.get_sample_name(d) for d in items]))
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
        if prioritize_by:
            work_dir = _sv_workdir(data)
            priority_files = [_prioritize_vcf(vcaller, vfile, prioritize_by, post_prior_fn, work_dir, data)
                              for vcaller, vfile, post_prior_fn in inputs]
            priority_tsv = _combine_files([xs[0] for xs in priority_files], work_dir, data)
            raw_files = {}
            for svcaller, fname in zip([xs[0] for xs in inputs], [xs[1] for xs in priority_files]):
                clean_fname = os.path.join(os.path.dirname(fname), "%s-%s-prioritized%s" %
                                           (dd.get_sample_name(data), svcaller, utils.splitext_plus(fname)[-1]))
                utils.symlink_plus(fname, clean_fname)
                raw_files[svcaller] = clean_fname
            data["sv"].append({"variantcaller": "sv-prioritize", "vrn_file": priority_tsv,
                               "raw_files": raw_files})
    # Disabled on move to CWL, not used and tested with CNVkit changes
    # data = _cnv_prioritize(data)
    return [data]

def is_gene_list(bed_file):
    """Check if the file is only a list of genes, not a BED
    """
    with utils.open_gzipsafe(bed_file) as in_handle:
        for line in in_handle:
            if not line.startswith("#"):
                if len(line.split()) == 1:
                    return True
                else:
                    return False

def _find_gene_list_from_bed(bed_file, base_file, data):
    """Retrieve list of gene names from input BED file.
    """
    # Check for a gene list, we can just return that.
    if is_gene_list(bed_file):
        return bed_file
    out_file = "%s-genes.txt" % utils.splitext_plus(base_file)[0]
    if not os.path.exists(out_file):
        genes = set([])
        import pybedtools
        with utils.open_gzipsafe(bed_file) as in_handle:
            for r in pybedtools.BedTool(in_handle):
                if r.name:
                    if not r.name.startswith("{"):
                        genes.add(r.name)
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                if len(genes) > 0:
                    out_handle.write("\n".join(sorted(list(genes))) + "\n")
    if utils.file_exists(out_file):
        return out_file

def _prioritize_vcf(caller, vcf_file, prioritize_by, post_prior_fn, work_dir, data):
    """Provide prioritized tab delimited output for a single caller.
    """
    sample = dd.get_sample_name(data)
    out_file = os.path.join(work_dir, "%s-%s-prioritize.tsv" % (sample, caller))
    simple_vcf = os.path.join(work_dir, "%s-%s-simple.vcf.gz" % (sample, caller))
    if not utils.file_exists(simple_vcf):
        gene_list = _find_gene_list_from_bed(prioritize_by, out_file, data)
        # If we have a standard gene list we can skip BED based prioritization
        priority_vcf = "%s.vcf.gz" % utils.splitext_plus(out_file)[0]
        if gene_list:
            if vcf_file.endswith(".vcf.gz"):
                utils.symlink_plus(vcf_file, priority_vcf)
            else:
                assert vcf_file.endswith(".vcf")
                utils.symlink_plus(vcf_file, priority_vcf.replace(".vcf.gz", ".vcf"))
                vcfutils.bgzip_and_index(priority_vcf.replace(".vcf.gz", ".vcf"),
                                         data["config"], remove_orig=False)
        # otherwise prioritize based on BED and proceed
        else:
            if not utils.file_exists(priority_vcf):
                with file_transaction(data, priority_vcf) as tx_out_file:
                    resources = config_utils.get_resources("bcbio_prioritize", data["config"])
                    jvm_opts = resources.get("jvm_opts", ["-Xms1g", "-Xmx4g"])
                    jvm_opts = config_utils.adjust_opts(jvm_opts, {"algorithm": {"memory_adjust":
                                                                                 {"direction": "increase",
                                                                                  "maximum": "30000M",
                                                                                  "magnitude": dd.get_cores(data)}}})
                    jvm_opts = " ".join(jvm_opts)
                    export = utils.local_path_export()
                    cmd = ("{export} bcbio-prioritize {jvm_opts} known -i {vcf_file} -o {tx_out_file} "
                           " -k {prioritize_by}")
                    do.run(cmd.format(**locals()), "Prioritize: select in known regions of interest")

        data_dir = os.path.dirname(os.path.realpath(utils.which("simple_sv_annotation.py")))
        with file_transaction(data, simple_vcf) as tx_out_file:
            fusion_file = os.path.join(data_dir, "fusion_pairs.txt")
            opts = ""
            if os.path.exists(fusion_file):
                opts += " --known_fusion_pairs %s" % fusion_file
            if not gene_list:
                opts += " --gene_list %s" % os.path.join(data_dir, "az-cancer-panel.txt")
            else:
                opts += " --gene_list %s" % gene_list
            cmd = "simple_sv_annotation.py {opts} -o - {priority_vcf} | bgzip -c > {tx_out_file}"
            do.run(cmd.format(**locals()), "Prioritize: simplified annotation output")
    simple_vcf = vcfutils.bgzip_and_index(vcfutils.sort_by_ref(simple_vcf, data), data["config"])
    if post_prior_fn:
        simple_vcf = post_prior_fn(simple_vcf, work_dir, data)
    if not utils.file_uptodate(out_file, simple_vcf):
        with file_transaction(data, out_file) as tx_out_file:
            export = utils.local_path_export(env_cmd="vawk")
            cmd = ("{export} zcat {simple_vcf} | vawk -v SNAME={sample} -v CALLER={caller} "
                   """'{{if (($7 == "PASS" || $7 == ".") && (S${sample}$GT != "0/0")) """
                   "print CALLER,SNAME,$1,$2,I$END,"
                   """I$SVTYPE=="BND" ? I$SVTYPE":"$3":"I$MATEID : I$SVTYPE,"""
                   "I$LOF,I$SIMPLE_ANN,"
                   "S${sample}$SR,S${sample}$PE,S${sample}$PR}}' > {tx_out_file}")
            do.run(cmd.format(**locals()), "Prioritize: convert to tab delimited")
    return out_file, simple_vcf

def _combine_files(tsv_files, work_dir, data):
    """Combine multiple priority tsv files into a final sorted output.
    """
    header = "\t".join(["caller", "sample", "chrom", "start", "end", "svtype",
                        "lof", "annotation", "split_read_support", "paired_support_PE", "paired_support_PR"])
    sample = dd.get_sample_name(data)
    out_file = os.path.join(work_dir, "%s-prioritize.tsv" % (sample))
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            tmpdir = os.path.dirname(tx_out_file)
            input_files = " ".join(tsv_files)
            sort_cmd = bedutils.get_sort_cmd(tmpdir)
            cmd = "{{ echo '{header}'; cat {input_files} | {sort_cmd} -k3,3 -k4,4n; }} > {tx_out_file}"
            do.run(cmd.format(**locals()), "Combine prioritized from multiple callers")
    return out_file

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "prioritize"))

# ## CNV prioritization by genes of interest and confidence intervals

def _cnvkit_prioritize(sample, genes, allele_file, metrics_file):
    """Summarize non-diploid calls with copy numbers and confidence intervals.
    """
    mdf = pd.read_csv(metrics_file, sep="\t")
    mdf.columns = [x.lower() for x in mdf.columns]
    if len(genes) > 0:
        mdf = mdf[mdf["gene"].str.contains("|".join(genes))]
    mdf = mdf[["chromosome", "start", "end", "gene", "log2", "ci_hi", "ci_lo"]]
    adf = pd.read_csv(allele_file, sep="\t")
    if len(genes) > 0:
        adf = adf[adf["gene"].str.contains("|".join(genes))]
    if "cn1" in adf.columns and "cn2" in adf.columns:
        adf = adf[["chromosome", "start", "end", "cn", "cn1", "cn2"]]
    else:
        adf = adf[["chromosome", "start", "end", "cn"]]
    df = pd.merge(mdf, adf, on=["chromosome", "start", "end"])
    df = df[df["cn"] != 2]
    if len(df) > 0:
        def passes(row):
            spread = abs(row["ci_hi"] - row["ci_lo"])
            return spread < 0.25
        df["passes"] = df.apply(passes, axis=1)
    df.insert(0, "sample", [sample] * len(df))
    return df

def _cnv_prioritize(data):
    """Perform confidence interval based prioritization for CNVs.
    """
    supported = {"cnvkit": {"inputs": ["call_file", "segmetrics"], "fn": _cnvkit_prioritize}}
    pcall = None
    priority_files = None
    for call in data.get("sv", []):
        if call["variantcaller"] in supported:
            priority_files = [call.get(x) for x in supported[call["variantcaller"]]["inputs"]]
            priority_files = [x for x in priority_files if x is not None and utils.file_exists(x)]
            if len(priority_files) == len(supported[call["variantcaller"]]["inputs"]):
                pcall = call
                break
    prioritize_by = tz.get_in(["config", "algorithm", "svprioritize"], data)
    if pcall and prioritize_by:
        out_file = "%s-prioritize.tsv" % utils.splitext_plus(priority_files[0])[0]
        gene_list = _find_gene_list_from_bed(prioritize_by, out_file, data)
        if gene_list:
            with open(gene_list) as in_handle:
                genes = [x.strip() for x in in_handle]
            args = [dd.get_sample_name(data), genes] + priority_files
            df = supported[pcall["variantcaller"]]["fn"](*args)
            with file_transaction(data, out_file) as tx_out_file:
                df.to_csv(tx_out_file, sep="\t", index=False)
            pcall["priority"] = out_file
    return data
