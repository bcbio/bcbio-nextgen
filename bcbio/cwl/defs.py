"""Definitions of workflows for translation into common workflow language.

This organizes the metadata and other information about workflows,
providing the necessary information to translate into CWL. The goal is to
eventually replace pipeline/main.py workflows with generalized
versions of this code.

Ideally we could translate the specification of these workflows into a YAML-based
DSL (some sort of CWL-lite), instead of writing them in Python code.

The variable 'workflows' provides a dictionary to retrieve the steps and outputs
for each of the defined workflows.
"""
import collections
from bcbio.pipeline import datadict as dd

def s(name, parallel, inputs, outputs, image, programs=None, disk=None, cores=None, unlist=None,
      no_files=False):
    """Represent a step in a workflow.

    name -- The run function name, which must match a definition in distributed/multitasks
    inputs -- List of input keys required for the function. Each key is of the type:
      ["toplevel", "sublevel"] -- an argument you could pass to toolz.get_in.
    outputs -- List of outputs with information about file type. Use cwlout functions
    programs -- Required programs for this step, used to define resource usage.
    disk -- Information about disk usage requirements, specified as multipliers of
            input files. Ensures enough disk present when that is a limiting factor
            when selecting cloud node resources.
    cores -- Maximum cores necessary for this step, for non-multicore processes.
    unlist -- Variables being unlisted by this process. Useful for parallelization splitting and
      batching from multiple variables, like variant calling.
    no_files -- This step does not require file access.
    parallel -- Parallelization approach. There are three different levels of parallelization,
      each with subcomponents:

      1. multi -- Multiple samples, parallelizing at the sample level. Used in top-level workflow.
        - multi-parallel -- Run individual samples in parallel.
        - multi-combined -- Run all samples together.
        - multi-batch -- Run all samples together, converting into batches of grouped samples.
      2. single -- A single sample, used in sub-workflows.
        - single-split -- Split a sample into sub-components (by read sections).
        - single-parallel -- Run sub-components of a sample in parallel.
        - single-merge -- Merge multiple sub-components into a single sample.
        - single-single -- Single sample, single item, nothing fancy.
      3. batch -- Several related samples (tumor/normal, or populations). Used in sub-workflows.
        - batch-split -- Split a batch of samples into sub-components (by genomic region).
        - batch-parallel -- Run sub-components of a batch in parallel.
        - batch-merge -- Merge sub-components back into a single batch.
        - batch-single -- Run on a single batch.
    """
    Step = collections.namedtuple("Step", "name parallel inputs outputs image programs disk cores unlist no_files")
    if programs is None: programs = []
    if unlist is None: unlist = []
    return Step(name, parallel, inputs, outputs, image, programs, disk, cores, unlist, no_files)

def w(name, parallel, workflow, internal):
    """A workflow, allowing specification of sub-workflows for nested parallelization.

    name and parallel are documented under the Step (s) function.
    workflow -- a list of Step tuples defining the sub-workflow
    internal -- variables used in the sub-workflow but not exposed to subsequent steps
    """
    Workflow = collections.namedtuple("Workflow", "name parallel workflow internal")
    return Workflow(name, parallel, workflow, internal)

def et(name, parallel, inputs, outputs, expression):
    """Represent an ExpressionTool that reorders inputs using javascript.
    """
    ExpressionTool = collections.namedtuple("ExpressionTool", "name inputs outputs expression parallel")
    return ExpressionTool(name, inputs, outputs, expression, parallel)

def cwlout(key, valtype=None, extensions=None, fields=None, exclude=None):
    """Definition of an output variable, defining the type and associated secondary files.
    """
    out = {"id": key}
    if valtype:
        out["type"] = valtype
    if fields:
        out["fields"] = fields
    if extensions:
        out["secondaryFiles"] = extensions
    if exclude:
        out["exclude"] = exclude
    return out

def _alignment(checkpoints):
    process_alignment_out = [cwlout(["work_bam"], ["File", "null"], [".bai"]),
                             cwlout(["align_bam"], ["File", "null"], [".bai"]),
                             cwlout(["hla", "fastq"], ["null", {"type": "array", "items": "File"}]),
                             cwlout(["work_bam_plus", "disc"], ["File", "null"], [".bai"]),
                             cwlout(["work_bam_plus", "sr"], ["File", "null"], [".bai"])]
    if checkpoints["umi"]:
        process_alignment_out.append(cwlout(["config", "algorithm", "rawumi_avg_cov"], ["int", "null"]))
        process_alignment_out.append(cwlout(["umi_bam"], ["File", "null"], [".bai"]))
    align = [s("prep_align_inputs", "single-split" if checkpoints["align_split"] else "single-single",
               [["alignment_rec"]],
               [cwlout("process_alignment_rec", "record",
                       fields=[cwlout(["files"], ["null", {"type": "array", "items": "File"}], [".gbi"]),
                               cwlout(["config", "algorithm", "quality_format"], ["string", "null"]),
                               cwlout(["align_split"], ["string", "null"])])],
               "bcbio-vc", ["grabix", "htslib", "biobambam", "atropos;env=python3",
                            "optitype;env=python2", "razers3=3.5.0", "coincbc"],  # HLA deps for general docker inclusion
               disk={"files": 1.5}),
             s("process_alignment", "single-parallel" if checkpoints["align_split"] else "single-single",
               [["alignment_rec"], ["process_alignment_rec"]], process_alignment_out,
               "bcbio-vc", ["bwa", "bwakit", "grabix", "minimap2", "novoalign", "snap-aligner=1.0dev.97",
                            "sentieon;env=python2", "samtools", "pysam>=0.13.0", "sambamba", "fgbio",
                            "umis;env=python2",
                            "biobambam", "seqtk", "samblaster", "variantbam"],
               disk={"files": 2})]
    if checkpoints["align_split"]:
        inp = [["alignment_rec"], ["work_bam"], ["align_bam"],
               ["work_bam_plus", "disc"], ["work_bam_plus", "sr"], ["hla", "fastq"]]
        outp = [cwlout(["align_bam"], ["File", "null"], [".bai"]),
                cwlout(["work_bam_plus", "disc"], ["File", "null"], [".bai"]),
                cwlout(["work_bam_plus", "sr"], ["File", "null"], [".bai"]),
                cwlout(["hla", "fastq"], ["null", {"type": "array", "items": "File"}])]
        if checkpoints["umi"]:
            inp.append(["umi_bam"])
            inp.append(["config", "algorithm", "rawumi_avg_cov"])
            outp.append(cwlout(["umi_bam"], ["File", "null"], [".bai"]))
            outp.append(cwlout(["config", "algorithm", "rawumi_avg_cov"], ["int", "null"]))
        align += [s("merge_split_alignments", "single-merge", inp, outp,
                  "bcbio-vc", ["biobambam", "samtools", "variantbam"],
                    disk={"files": 3.5})]
    return align

def _variant_hla(checkpoints):
    """Add hla analysis to workflow, if configured.
    """
    if not checkpoints.get("hla"):
        return [], []
    hla = [s("hla_to_rec", "multi-batch",
             [["hla", "fastq"],
              ["config", "algorithm", "hlacaller"]],
               [cwlout("hla_rec", "record")],
               "bcbio-vc", cores=1, no_files=True),
           s("call_hla", "multi-parallel",
             [["hla_rec"]],
             [cwlout(["hla", "hlacaller"], ["string", "null"]),
              cwlout(["hla", "call_file"], ["File", "null"])],
             "bcbio-vc", ["optitype;env=python2", "razers3=3.5.0", "coincbc"])]
    return hla, [["hla", "call_file"]]

def _variant_vc(checkpoints):
    """Add variant calling to workflow, if configured.
    """
    if not checkpoints.get("vc"):
        return [], []
    vc_wf = [s("get_parallel_regions", "batch-split",
               [["batch_rec"]],
               [cwlout(["region_block"], {"type": "array", "items": "string"})],
               "bcbio-vc",
               disk={"files": 2.0}, cores=1),
             s("variantcall_batch_region", "batch-parallel",
               [["batch_rec"], ["region_block"]],
               [cwlout(["vrn_file_region"], ["File", "null"], [".tbi"]),
                cwlout(["region_block"], {"type": "array", "items": "string"})],
               "bcbio-vc", ["bcftools", "bedtools", "freebayes=1.1.0.46",
                            "gatk4", "vqsr_cnn", "deepvariant;env=dv", "sentieon;env=python2",
                            "htslib", "octopus", "picard", "platypus-variant;env=python2", "pythonpy",
                            "samtools", "pysam>=0.13.0", "strelka;env=python2", "vardict", "vardict-java",
                            "varscan", "moreutils", "vcfanno", "vcflib", "vt", "r=3.5.1", "r-base",
                            "perl"],
               disk={"files": 2.0}),
             s("concat_batch_variantcalls", "batch-merge",
               [["batch_rec"], ["region_block"], ["vrn_file_region"]],
               [cwlout(["vrn_file"], "File", [".tbi"])],
               "bcbio-vc", ["bcftools", "htslib", "gatk4"],
               disk={"files": 1.5}, cores=1)]
    if not checkpoints.get("jointvc"):
        vc_wf += [s("postprocess_variants", "batch-single",
                    [["batch_rec"], ["vrn_file"]],
                    [cwlout(["vrn_file"], "File", [".tbi"])],
                    "bcbio-vc", ["snpeff=4.3.1t"], disk={"files": 0.5})]
    vc_rec_exclude = [["align_bam"]]
    if not checkpoints.get("jointvc"):
        vc_rec_exclude.append(["genome_resources", "variation"])
    vc_wf += [s("compare_to_rm", "batch-single",
                [["batch_rec"], ["vrn_file"]],
                [cwlout("vc_rec", "record",
                        fields=[cwlout(["batch_samples"], ["null", {"type": "array", "items": "string"}]),
                                cwlout(["validate", "summary"], ["File", "null"]),
                                cwlout(["validate", "tp"], ["File", "null"], [".tbi"]),
                                cwlout(["validate", "fp"], ["File", "null"], [".tbi"]),
                                cwlout(["validate", "fn"], ["File", "null"], [".tbi"]),
                                cwlout("inherit", exclude=vc_rec_exclude)])],
                "bcbio-vc", ["bcftools", "bedtools", "pythonpy", "gvcf-regions;env=python2",
                             "htslib", "rtg-tools", "vcfanno"],
                disk={"files": 1.5})]
    batch_in = [["analysis"], ["genome_build"], ["align_bam"], ["vrn_file"],
                ["metadata", "batch"], ["metadata", "phenotype"],
                ["config", "algorithm", "callable_regions"], ["regions", "sample_callable"],
                ["config", "algorithm", "variantcaller"],
                ["config", "algorithm", "ensemble"],
                ["config", "algorithm", "vcfanno"],
                ["config", "algorithm", "coverage_interval"],
                ["config", "algorithm", "effects"],
                ["config", "algorithm", "min_allele_fraction"],
                ["config", "algorithm", "exclude_regions"],
                ["config", "algorithm", "variant_regions"],
                ["config", "algorithm", "variant_regions_merged"],
                ["config", "algorithm", "validate"], ["config", "algorithm", "validate_regions"],
                ["config", "algorithm", "tools_on"],
                ["config", "algorithm", "tools_off"],
                ["reference", "fasta", "base"],
                ["reference", "rtg"], ["reference", "genome_context"],
                ["genome_resources", "variation", "clinvar"],
                ["genome_resources", "variation", "cosmic"], ["genome_resources", "variation", "dbsnp"],
                ["genome_resources", "variation", "esp"], ["genome_resources", "variation", "exac"],
                ["genome_resources", "variation", "gnomad_exome"],
                ["genome_resources", "variation", "1000g"],
                ["genome_resources", "variation", "lcr"], ["genome_resources", "variation", "polyx"],
                ["genome_resources", "variation", "encode_blacklist"],
                ["genome_resources", "aliases", "ensembl"], ["genome_resources", "aliases", "human"],
                ["genome_resources", "aliases", "snpeff"], ["reference", "snpeff", "genome_build"]]
    if checkpoints.get("umi"):
        batch_in.append(["config", "algorithm", "umi_type"])
    if checkpoints.get("rnaseq"):
        batch_in += [["genome_resources", "variation", "editing"]]
    else:
        batch_in += [["genome_resources", "variation", "train_hapmap"],
                     ["genome_resources", "variation", "train_indels"]]
    vc = [s("batch_for_variantcall", "multi-batch", batch_in,
            [cwlout("batch_rec", "record",
                    fields=[cwlout(["config", "algorithm", "variantcaller_order"], "int"),
                            cwlout("inherit")])],
            "bcbio-vc",
            disk={"files": 2.0}, cores=1,
            unlist=[["config", "algorithm", "variantcaller"]], no_files=True),
          w("variantcall", "multi-parallel", vc_wf,
            [["region"], ["region_block"], ["vrn_file_region"], ["vrn_file"], ["validate", "summary"]])]
    if checkpoints.get("jointvc"):
        vc += _variant_jointvc()
    if checkpoints.get("ensemble"):
        vc += _variant_ensemble(checkpoints)
    summarize_in = [["jointvc_rec" if checkpoints.get("jointvc") else "vc_rec"]]
    if checkpoints.get("ensemble"):
        summarize_in += [["ensemble_rec"]]
    vc += [s("summarize_vc", "multi-combined", summarize_in,
             [cwlout(["variants", "calls"], {"type": "array", "items": ["File", "null"]}),
              cwlout(["variants", "gvcf"], ["null", {"type": "array", "items": ["File", "null"]}]),
              cwlout(["variants", "samples"], {"type": "array", "items": {"type": "array",
                                                                          "items": ["File", "null"]}}),
              cwlout(["validate", "grading_summary"], ["File", "null"]),
              cwlout(["validate", "grading_plots"], {"type": "array", "items": ["File", "null"]})],
             "bcbio-vc",
             disk={"files": 2.0}, cores=1)]
    return vc, [["validate", "grading_summary"], ["variants", "calls"], ["variants", "gvcf"]]

def _variant_ensemble(checkpoints):
    out = [s("batch_for_ensemble", "multi-combined",
             [["jointvc_rec" if checkpoints.get("jointvc") else "vc_rec"]],
             [cwlout("ensemble_prep_rec", "record",
                     fields=[cwlout(["batch_id"], "string"),
                             cwlout(["variants", "calls"], {"type": "array", "items": "File"}),
                             cwlout(["variants", "variantcallers"], {"type": "array", "items": "string"}),
                             cwlout("inherit")])],
             "bcbio-vc", cores=1, no_files=True),
           s("combine_calls", "multi-parallel",
             ["ensemble_prep_rec"],
             [cwlout("ensemble_rec", "record",
                     fields=[cwlout(["ensemble", "vrn_file"], ["File", "null"]),
                             cwlout(["ensemble", "validate", "summary"], ["File", "null"]),
                             cwlout(["ensemble", "batch_samples"], {"type": "array", "items": "string"}),
                             cwlout(["ensemble", "batch_id"], "string")])],
             "bcbio-vc", ["bcbio-variation-recall"])]
    return out

def _variant_jointvc():
    wf = [s("get_parallel_regions_jointvc", "batch-split",
            [["jointvc_batch_rec"]],
            [cwlout(["region"], "string")],
            "bcbio-vc",
            disk={"files": 1.5}, cores=1),
          s("run_jointvc", "batch-parallel",
            [["jointvc_batch_rec"], ["region"]],
            [cwlout(["vrn_file_region"], ["File", "null"], [".tbi"]), cwlout(["region"], "string")],
            "bcbio-vc", ["gatk4", "gvcfgenotyper", "sentieon;env=python2"],
            disk={"files": 1.5}, cores=1),
          s("concat_batch_variantcalls_jointvc", "batch-merge",
            [["jointvc_batch_rec"], ["region"], ["vrn_file_region"]],
            [cwlout(["vrn_file_joint"], "File", [".tbi"])],
            "bcbio-vc", ["bcftools", "htslib", "gatk4"],
            disk={"files": 1.5}, cores=1),
          s("postprocess_variants", "batch-single",
            [["jointvc_batch_rec"], ["vrn_file_joint"]],
            [cwlout(["vrn_file_joint"], "File", [".tbi"])],
            "bcbio-vc", ["snpeff=4.3.1t"],
            disk={"files": 1.5}),
          s("finalize_jointvc", "batch-single",
            [["jointvc_batch_rec"], ["vrn_file_joint"]],
            [cwlout("jointvc_rec", "record")],
            "bcbio-vc",
            disk={"files": 1.5}, cores=1)]
    out = [s("batch_for_jointvc", "multi-batch",
             ["vc_rec"],
             [cwlout("jointvc_batch_rec", "record")],
             "bcbio-vc",
             disk={"files": 1.5}, cores=1, no_files=True),
           w("jointcall", "multi-parallel", wf,
             [["region"], ["vrn_file_region"], ["vrn_file"]])]
    return out

def _variant_checkpoints(samples):
    """Check sample configuration to identify required steps in analysis.
    """
    checkpoints = {}
    checkpoints["vc"] = any([dd.get_variantcaller(d) or d.get("vrn_file") for d in samples])
    checkpoints["sv"] = any([dd.get_svcaller(d) for d in samples])
    checkpoints["jointvc"] = any([(dd.get_jointcaller(d) or "gvcf" in dd.get_tools_on(d))
                                  for d in samples])
    checkpoints["hla"] = any([dd.get_hlacaller(d) for d in samples])
    checkpoints["align"] = any([(dd.get_aligner(d) or dd.get_bam_clean(d)) for d in samples])
    checkpoints["align_split"] = not all([(dd.get_align_split_size(d) is False or
                                           not dd.get_aligner(d))
                                          for d in samples])
    checkpoints["archive"] = any([dd.get_archive(d) for d in samples])
    checkpoints["umi"] = any([dd.get_umi_consensus(d) for d in samples])
    checkpoints["ensemble"] = any([dd.get_ensemble(d) for d in samples])
    checkpoints["cancer"] = any(dd.get_phenotype(d) in ["tumor"] for d in samples)
    return checkpoints

def _postprocess_alignment(checkpoints):
    wf = [s("prep_samples_to_rec", "multi-combined",
            [["config", "algorithm", "coverage"],
             ["rgnames", "sample"],
             ["config", "algorithm", "background", "cnv_reference"],
             ["config", "algorithm", "svcaller"],
             ["config", "algorithm", "sv_regions"],
             ["config", "algorithm", "variant_regions"],
             ["reference", "fasta", "base"]],
            [cwlout("prep_samples_rec", "record")],
            "bcbio-vc",
            disk={"files": 0.5}, cores=1, no_files=True),
          s("prep_samples", "multi-parallel",
            ["prep_samples_rec"],
            [cwlout(["rgnames", "sample"], "string"),
             cwlout(["config", "algorithm", "variant_regions"], ["File", "null"]),
             cwlout(["config", "algorithm", "variant_regions_merged"], ["File", "null"]),
             cwlout(["config", "algorithm", "variant_regions_orig"], ["File", "null"]),
             cwlout(["config", "algorithm", "coverage"], ["File", "null"]),
             cwlout(["config", "algorithm", "coverage_merged"], ["File", "null"]),
             cwlout(["config", "algorithm", "coverage_orig"], ["File", "null"]),
             cwlout(["config", "algorithm", "seq2c_bed_ready"], ["File", "null"])],
            "bcbio-vc", ["htslib", "bedtools", "pythonpy"],
            disk={"files": 0.5}, cores=1),
          s("postprocess_alignment_to_rec", "multi-combined",
            [["align_bam"],
             ["config", "algorithm", "archive"],
             ["config", "algorithm", "coverage_interval"],
             ["config", "algorithm", "exclude_regions"],
             ["config", "algorithm", "variant_regions"],
             ["config", "algorithm", "variant_regions_merged"],
             ["config", "algorithm", "variant_regions_orig"],
             ["config", "algorithm", "coverage"],
             ["config", "algorithm", "coverage_merged"],
             ["config", "algorithm", "coverage_orig"],
             ["config", "algorithm", "seq2c_bed_ready"],
             ["config", "algorithm", "recalibrate"],
             ["config", "algorithm", "tools_on"],
             ["genome_resources", "rnaseq", "gene_bed"],
             ["genome_resources", "variation", "dbsnp"],
             ["genome_resources", "variation", "lcr"], ["genome_resources", "variation", "polyx"],
             ["genome_resources", "variation", "encode_blacklist"],
             ["reference", "fasta", "base"]],
            [cwlout("postprocess_alignment_rec", "record")],
            "bcbio-vc",
            disk={"files": 1.5}, cores=1, no_files=True),
          s("postprocess_alignment", "multi-parallel",
            [["postprocess_alignment_rec"]],
            [cwlout(["config", "algorithm", "coverage_interval"], ["string", "null"]),
             cwlout(["config", "algorithm", "variant_regions"], ["File", "null"]),
             cwlout(["config", "algorithm", "variant_regions_merged"], ["File", "null"]),
             cwlout(["config", "algorithm", "variant_regions_orig"], ["File", "null"]),
             cwlout(["config", "algorithm", "coverage"], ["File", "null"]),
             cwlout(["config", "algorithm", "coverage_merged"], ["File", "null"]),
             cwlout(["config", "algorithm", "coverage_orig"], ["File", "null"]),
             cwlout(["config", "algorithm", "seq2c_bed_ready"], ["File", "null"]),
             cwlout(["regions", "callable"], ["File", "null"]),
             cwlout(["regions", "sample_callable"], ["File", "null"]),
             cwlout(["regions", "nblock"], ["File", "null"]),
             cwlout(["depth", "samtools", "stats"], ["File", "null"]),
             cwlout(["depth", "samtools", "idxstats"], ["File", "null"]),
             cwlout(["depth", "variant_regions", "regions"], ["File", "null"]),
             cwlout(["depth", "variant_regions", "dist"], ["File", "null"]),
             cwlout(["depth", "sv_regions", "regions"], ["File", "null"]),
             cwlout(["depth", "sv_regions", "dist"], ["File", "null"]),
             cwlout(["depth", "coverage", "regions"], ["File", "null"]),
             cwlout(["depth", "coverage", "dist"], ["File", "null"]),
             cwlout(["depth", "coverage", "thresholds"], ["File", "null"]),
             cwlout(["align_bam"], ["File", "null"])],
            "bcbio-vc", ["sambamba", "goleft", "bedtools", "htslib", "gatk4", "mosdepth", "sentieon;env=python2"],
            disk={"files": 3.0}),
          s("combine_sample_regions", "multi-combined",
            [["regions", "callable"], ["regions", "nblock"], ["metadata", "batch"],
             ["config", "algorithm", "nomap_split_size"], ["config", "algorithm", "nomap_split_targets"],
             ["reference", "fasta", "base"]],
            [cwlout(["config", "algorithm", "callable_regions"], "File"),
             cwlout(["config", "algorithm", "non_callable_regions"], "File"),
             cwlout(["config", "algorithm", "callable_count"], "int")],
            "bcbio-vc", ["bedtools", "htslib", "gatk4"],
            disk={"files": 0.5}, cores=1)]
    out = [["regions", "sample_callable"]]
    if checkpoints.get("archive"):
        wf += [s("archive_to_cram", "multi-parallel",
                 [["postprocess_alignment_rec"]],
                 [cwlout(["archive_bam"], ["File", "null"], [".crai"])],
                 "bcbio-vc", ["samtools"],
                 disk={"files": 3.0})]
        out += [["archive_bam"]]
    return wf, out

def variant(samples):
    """Variant calling workflow definition for CWL generation.
    """
    checkpoints = _variant_checkpoints(samples)
    if checkpoints["align"]:
        align_wf = _alignment(checkpoints)
        alignin = [["files"], ["analysis"],
                   ["config", "algorithm", "align_split_size"],
                   ["reference", "fasta", "base"],
                   ["rgnames", "pl"], ["rgnames", "sample"], ["rgnames", "pu"],
                   ["rgnames", "lane"], ["rgnames", "rg"], ["rgnames", "lb"],
                   ["reference", "aligner", "indexes"],
                   ["config", "algorithm", "aligner"],
                   ["config", "algorithm", "trim_reads"],
                   ["config", "algorithm", "adapters"],
                   ["config", "algorithm", "bam_clean"],
                   ["config", "algorithm", "variant_regions"],
                   ["config", "algorithm", "mark_duplicates"]]
        if checkpoints["hla"]:
            alignin.append(["config", "algorithm", "hlacaller"])
        if checkpoints["umi"]:
            alignin.append(["config", "algorithm", "umi_type"])
        align = [s("alignment_to_rec", "multi-combined", alignin,
                   [cwlout("alignment_rec", "record")],
                   "bcbio-vc",
                   disk={"files": 1.5}, cores=1, no_files=True),
                 w("alignment", "multi-parallel", align_wf,
                   [["align_split"], ["process_alignment_rec"],
                    ["work_bam"], ["config", "algorithm", "quality_format"]])]
    else:
        align = [s("organize_noalign", "multi-parallel",
                   ["files"],
                   [cwlout(["align_bam"], ["File", "null"], [".bai"]),
                    cwlout(["work_bam_plus", "disc"], ["File", "null"]),
                    cwlout(["work_bam_plus", "sr"], ["File", "null"]),
                    cwlout(["hla", "fastq"], ["File", "null"])],
                   "bcbio-vc", cores=1)]
    align_out = [["rgnames", "sample"], ["align_bam"]]
    pp_align, pp_align_out = _postprocess_alignment(checkpoints)
    if checkpoints["umi"]:
        align_out += [["umi_bam"]]
    vc, vc_out = _variant_vc(checkpoints)
    sv, sv_out = _variant_sv(checkpoints)
    hla, hla_out = _variant_hla(checkpoints)
    qc, qc_out = _qc_workflow(checkpoints)
    steps = align + pp_align + hla + vc + sv + qc
    final_outputs = align_out + pp_align_out + vc_out + hla_out + sv_out + qc_out
    return steps, final_outputs

def _qc_workflow(checkpoints):
    qc_inputs = \
      [["align_bam"], ["analysis"], ["reference", "fasta", "base"], ["reference", "versions"],
       ["config", "algorithm", "tools_on"], ["config", "algorithm", "tools_off"],
       ["genome_build"], ["config", "algorithm", "qc"], ["metadata", "batch"], ["metadata", "phenotype"],
       ["config", "algorithm", "coverage_interval"],
       ["depth", "variant_regions", "regions"], ["depth", "variant_regions", "dist"],
       ["depth", "samtools", "stats"], ["depth", "samtools", "idxstats"],
       ["depth", "sv_regions", "regions"], ["depth", "sv_regions", "dist"],
       ["depth", "coverage", "regions"], ["depth", "coverage", "dist"], ["depth", "coverage", "thresholds"],
       ["config", "algorithm", "variant_regions"],
       ["config", "algorithm", "variant_regions_merged"],
       ["config", "algorithm", "coverage"],
       ["config", "algorithm", "coverage_merged"]]
    if checkpoints.get("vc"):
        qc_inputs += [["variants", "samples"]]
    if checkpoints.get("umi"):
        qc_inputs += [["config", "algorithm", "umi_type"], ["config", "algorithm", "rawumi_avg_cov"], ["umi_bam"]]
    if checkpoints.get("cancer"):
        qc_inputs += [["reference", "viral"]]
    qc = [s("qc_to_rec", "multi-combined",
            qc_inputs, [cwlout("qc_rec", "record")],
            "bcbio-vc", disk={"files": 1.5}, cores=1, no_files=True),
          s("pipeline_summary", "multi-parallel",
            ["qc_rec"],
            [cwlout("qcout_rec", "record",
                    fields=[cwlout(["summary", "qc"], ["File", "null"]),
                            cwlout(["summary", "metrics"], ["string", "null"]),
                            cwlout(["genome_build"]), cwlout(["description"]),
                            cwlout(["reference", "versions"]),
                            cwlout(["config", "algorithm", "tools_off"]),
                            cwlout(["config", "algorithm", "tools_on"]),
                            cwlout(["config", "algorithm", "qc"])])],
            "bcbio-vc", ["bcftools", "bedtools", "fastqc=0.11.7=5", "goleft", "hts-nim-tools", "mosdepth",
                         "picard", "pythonpy", "qsignature", "qualimap", "sambamba",
                         "samtools", "preseq", "peddy;env=python2", "verifybamid2"],
            disk={"files": 2.0}),
          s("multiqc_summary", "multi-combined",
            [["qcout_rec"]],
            [cwlout(["summary", "multiqc"], ["File", "null"]),
             cwlout(["versions", "tools"], ["File", "null"]),
             cwlout(["versions", "data"], ["File", "null"])],
            "bcbio-vc", ["multiqc", "multiqc-bcbio"],
            disk={"files": 2.0}, cores=1)]
    qc_out = [cwlout(["summary", "multiqc"], {"type": "array", "items": ["File", "null"]}),
              ["versions", "tools"], ["versions", "data"]]
    return qc, qc_out

def _variant_sv(checkpoints):
    """Structural variant workflow.
    """
    if not checkpoints.get("sv"):
        return [], []
    sv = [s("detect_sv", "batch-single",
            [["sv_batch_rec"]],
            [cwlout("sv_rec", "record",
                    fields=[cwlout(["sv", "variantcaller"], ["string", "null"]),
                            cwlout(["sv", "vrn_file"], ["File", "null"], [".tbi"]),
                            cwlout(["sv", "supplemental"], {"type": "array", "items": ["File"]}),
                            cwlout(["svvalidate", "summary"], ["File", "null"]),
                            cwlout("inherit", exclude=[["align_bam"], ["work_bam_plus"],
                                                       ["reference", "snpeff"]])])],
            "bcbio-vc", ["bedtools", "cnvkit", "delly", "duphold", "extract-sv-reads", "gsort",
                         "lumpy-sv;env=python2", "manta;env=python2", "break-point-inspector", "mosdepth", "samtools",
                         "smoove;env=python2", "pysam>=0.13.0",
                         "seq2c", "simple_sv_annotation;env=python2", "survivor", "svtools;env=python2",
                         "svtyper;env=python2",
                         "r=3.5.1", "r-base", "xorg-libxt", "vawk;env=python2"],
            disk={"files": 2.0})]
    sv_batch_inputs = [["analysis"], ["genome_build"],
                       ["work_bam_plus", "disc"], ["work_bam_plus", "sr"],
                       ["config", "algorithm", "background", "cnv_reference"],
                       ["config", "algorithm", "tools_on"],
                       ["config", "algorithm", "tools_off"],
                       ["config", "algorithm", "svprioritize"],
                       ["config", "algorithm", "svvalidate"], ["regions", "sample_callable"],
                       ["genome_resources", "variation", "gc_profile"],
                       ["genome_resources", "variation", "germline_het_pon"],
                       ["genome_resources", "aliases", "snpeff"], ["reference", "snpeff", "genome_build"],
                       ["sv_coverage_rec"]]
    if checkpoints.get("vc"):
        sv_batch_inputs.append(["variants", "samples"])
    steps = [s("calculate_sv_bins", "multi-combined",
               [["align_bam"], ["reference", "fasta", "base"],
                ["metadata", "batch"], ["metadata", "phenotype"],
                ["config", "algorithm", "background", "cnv_reference"],
                ["config", "algorithm", "callable_regions"],
                ["config", "algorithm", "coverage_interval"],
                ["config", "algorithm", "exclude_regions"],
                ["config", "algorithm", "sv_regions"],
                ["config", "algorithm", "variant_regions"],
                ["config", "algorithm", "variant_regions_merged"],
                ["config", "algorithm", "seq2c_bed_ready"],
                ["config", "algorithm", "svcaller"],
                ["depth", "variant_regions", "regions"],
                ["genome_resources", "variation", "lcr"], ["genome_resources", "variation", "polyx"],
                ["genome_resources", "variation", "encode_blacklist"],
                ["genome_resources", "rnaseq", "gene_bed"]],
               [cwlout("sv_bin_rec", "record",
                       fields=[cwlout(["regions", "bins", "target"], ["File", "null"]),
                               cwlout(["regions", "bins", "antitarget"], ["File", "null"]),
                               cwlout(["regions", "bins", "gcannotated"], ["File", "null"]),
                               cwlout(["regions", "bins", "group"], ["string", "null"]),
                               cwlout("inherit")])],
               "bcbio-vc", ["bedtools", "cnvkit"],
               disk={"files": 1.5}, cores=1),
             s("calculate_sv_coverage", "multi-parallel",
               [["sv_bin_rec"]],
               [cwlout("sv_rawcoverage_rec", "record",
                       fields=[cwlout(["depth", "bins", "target"], ["File", "null"]),
                               cwlout(["depth", "bins", "antitarget"], ["File", "null"]),
                               cwlout(["depth", "bins", "seq2c"], ["File", "null"]),
                               cwlout("inherit")])],
               "bcbio-vc", ["mosdepth", "cnvkit", "seq2c"],
               disk={"files": 1.5}),
             s("normalize_sv_coverage", "multi-combined",
               [["sv_rawcoverage_rec"]],
               [cwlout("sv_coverage_rec", "record",
                       fields=[cwlout(["depth", "bins", "normalized"], ["File", "null"]),
                               cwlout(["depth", "bins", "background"], ["File", "null"]),
                               cwlout("inherit")])],
               "bcbio-vc", ["cnvkit"],
               disk={"files": 1.5}),
             s("batch_for_sv", "multi-batch", sv_batch_inputs,
               [cwlout("sv_batch_rec", "record")],
               "bcbio-vc",
               unlist=[["config", "algorithm", "svcaller"]]),
             w("svcall", "multi-parallel", sv, []),
             s("summarize_sv", "multi-combined",
               [["sv_rec"]],
               [cwlout(["sv", "calls"], {"type": "array", "items": ["File", "null"]}),
                cwlout(["sv", "supplemental"], {"type": "array", "items": ["File"]}),
                cwlout(["sv", "prioritize", "tsv"], {"type": "array", "items": ["File", "null"]}),
                cwlout(["sv", "prioritize", "raw"], {"type": "array", "items": ["File", "null"]}),
                cwlout(["svvalidate", "grading_summary"], ["File", "null"]),
                cwlout(["svvalidate", "grading_plots"], {"type": "array", "items": ["File", "null"]})],
               "bcbio-vc", ["bcbio-prioritize"], disk={"files": 1.0}, cores=1)]
    final_outputs = [["sv", "calls"], ["svvalidate", "grading_summary"], ["sv", "prioritize", "tsv"],
                     ["sv", "prioritize", "raw"], ["sv", "supplemental"]]
    return steps, final_outputs

def rnaseq(samples):
    checkpoints = _rnaseq_checkpoints(samples)
    prep = [s("prepare_sample", "multi-parallel",
              [["files"], dd.get_keys("sample_name"),
               dd.get_keys("ref_file"), dd.get_keys("genome_build"), dd.get_keys("gtf_file"),
               ["analysis"],
               ["rgnames", "pl"], ["rgnames", "pu"], ["rgnames", "lane"],
               ["rgnames", "rg"], ["rgnames", "lb"],
               ["reference", "aligner", "indexes"],
               ["config", "algorithm", "aligner"],
               ["config", "algorithm", "expression_caller"],
               ["config", "algorithm", "fusion_caller"],
               ["config", "algorithm", "quality_format"]],
              [cwlout("prep_rec", "record")],
              "bcbio-rnaseq", programs=["picard", "samtools", "pysam>=0.13.0"]),
            s("trim_sample", "multi-parallel",
              [["prep_rec"]],
              [cwlout("trim_rec", "record")],
              "bcbio-rnaseq", programs=["atropos;env=python3"])]
    align = [s("process_alignment", "multi-parallel",
               [["trim_rec"]],
               [cwlout(["align_bam"], "File", [".bai"])],
               "bcbio-rnaseq", ["star", "hisat2", "tophat;env=python2", "samtools",
                                "sambamba", "seqtk"],
               {"files": 1.5})]
    if checkpoints.get("vc"):
        pp_align, pp_align_out = _postprocess_alignment(checkpoints)
    else:
        pp_align, pp_align_out = [], []
    quantitate = [s("rnaseq_quantitate", "multi-parallel",
                  [["trim_rec"], ["align_bam"]],
                  [cwlout(dd.get_keys("count_file"), "File"),
                   cwlout(["quant", "tsv"], "File"),
                   cwlout(["quant", "hdf5"], "File"),
                   cwlout(["quant", "fusion"], "File")],
                  "bcbio-rnaseq", programs=["sailfish", "salmon", "kallisto>=0.43.1", "subread", "gffread",
                                            "r=3.5.1", "r-base", "xorg-libxt", "r-wasabi"],
                  disk={"files": 0.5})]
    qc = [s("qc_to_rec", "multi-combined",
            [["align_bam"], ["analysis"], ["reference", "fasta", "base"], dd.get_keys("gtf_file"),
             ["genome_build"], ["config", "algorithm", "coverage_interval"],
             ["config", "algorithm", "tools_on"], ["config", "algorithm", "tools_off"],
             ["config", "algorithm", "qc"]],
            [cwlout("qc_rec", "record")],
            "bcbio-rnaseq", disk={"files": 1.5}, cores=1, no_files=True),
          s("pipeline_summary", "multi-parallel",
            [["qc_rec"]],
            [cwlout("qcout_rec", "record",
                    fields=[cwlout(["summary", "qc"], ["File", "null"]),
                            cwlout(["summary", "metrics"], ["string", "null"]),
                            cwlout("inherit")])],
            "bcbio-rnaseq", ["bedtools", "fastqc=0.11.7=5", "goleft", "hts-nim-tools", "mosdepth",
                             "picard", "pythonpy", "qsignature", "qualimap",
                             "sambamba", "samtools"]),
          s("multiqc_summary", "multi-combined",
            [["qcout_rec"]],
            [cwlout(["summary", "multiqc"], ["File", "null"])],
            "bcbio-rnaseq", ["multiqc", "multiqc-bcbio"])]
    vc, vc_out = _variant_vc(checkpoints)
    fusion = [s("detect_fusions", "multi-parallel",
                [["quant", "fusion"], ["quant", "hdf5"], ["trim_rec"]],
                [cwlout(["fusion", "fasta"], "File"),
                 cwlout(["fusion", "json"], "File")],
                "bcbio-rnaseq", ["pizzly"])]

    steps = prep + align + pp_align + quantitate + qc + vc + fusion
    final_outputs = [["rgnames", "sample"], dd.get_keys("align_bam"), ["quant", "tsv"], ["summary", "multiqc"]] + \
                    vc_out + pp_align_out
    return steps, final_outputs

def _rnaseq_checkpoints(samples):
    """Check sample configuration to identify required steps in analysis.
    """
    checkpoints = {}
    checkpoints["rnaseq"] = True
    checkpoints["vc"] = any([dd.get_variantcaller(d) for d in samples])
    return checkpoints

workflows = \
  {"variant": variant, "variant2": variant, "rna-seq": rnaseq}
