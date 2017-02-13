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

def s(name, parallel, inputs, outputs, programs=None, disk=None, unlist=None):
    """Represent a step in a workflow.

    name -- The run function name, which must match a definition in distributed/multitasks
    inputs -- List of input keys required for the function. Each key is of the type:
      ["toplevel", "sublevel"] -- an argument you could pass to toolz.get_in.
    outputs -- List of outputs with information about file type. Use cwlout functions
    programs -- Required programs for this step, used to define resource usage.
    disk -- Information about disk usage requirements, specified as multipliers of
            input files. Ensures enough disk present when that is a limiting factor
            when selecting cloud node resources.
    unlist -- Variables being unlisted by this process. Useful for parallelization splitting and
      batching from multiple variables, like variant calling.
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
      3. batch -- Several related samples (tumor/normal, or populations). Used in sub-workflows.
        - batch-split -- Split a batch of samples into sub-components (by genomic region).
        - batch-parallel -- Run sub-components of a batch in parallel.
        - batch-merge -- Merge sub-components back into a single batch.
        - batch-single -- Run on a single batch.
    """
    Step = collections.namedtuple("Step", "name parallel inputs outputs programs disk unlist")
    if programs is None: programs = []
    if unlist is None: unlist = []
    return Step(name, parallel, inputs, outputs, programs, disk, unlist)

def w(name, parallel, workflow, internal):
    """A workflow, allowing specification of sub-workflows for nested parallelization.

    name and parallel are documented under the Step (s) function.
    workflow -- a list of Step tuples defining the sub-workflow
    internal -- variables used in the sub-workflow but not exposed to subsequent steps
    """
    Workflow = collections.namedtuple("Workflow", "name parallel workflow internal")
    return Workflow(name, parallel, workflow, internal)

def cwlout(key, valtype, extensions=None):
    """Definition of an output variable, defining the type and associated secondary files.
    """
    out = {"id": key,
           "type": valtype}
    if extensions:
        out["secondaryFiles"] = extensions
    return out

def _variant_shared():
    align = [s("prep_align_inputs", "single-split",
               [["files"],
                ["config", "algorithm", "align_split_size"],
                ["config", "algorithm", "aligner"]],
               [cwlout(["files"], "File", [".gbi"]),
                cwlout(["config", "algorithm", "quality_format"], "string"),
                cwlout(["align_split"], ["string", "null"])],
               ["bgzip", "pbgzip"],
               {"files": 1.5}),
             s("process_alignment", "single-parallel",
               [["files"], ["reference", "fasta", "base"], ["align_split"],
                ["rgnames", "pl"], ["rgnames", "sample"], ["rgnames", "pu"],
                ["rgnames", "lane"], ["rgnames", "rg"], ["rgnames", "lb"],
                ["reference", "aligner", "indexes"],
                ["config", "algorithm", "aligner"],
                ["config", "algorithm", "mark_duplicates"],
                ["config", "algorithm", "quality_format"]],
               [cwlout(["work_bam"], "File"),
                cwlout(["align_bam"], "File"),
                cwlout(["hla", "fastq"], ["File", "null"]),
                cwlout(["work_bam_plus", "disc"], ["File", "null"], [".bai"]),
                cwlout(["work_bam_plus", "sr"], ["File", "null"], [".bai"])],
               ["aligner", "samtools", "sambamba"],
               {"files": 1.5}),
             s("merge_split_alignments", "single-merge",
               [["work_bam"], ["align_bam"], ["work_bam_plus", "disc"], ["work_bam_plus", "sr"],
                ["hla", "fastq"]],
               [cwlout(["align_bam"], "File", [".bai"]),
                cwlout(["work_bam_plus", "disc"], ["File", "null"], [".bai"]),
                cwlout(["work_bam_plus", "sr"], ["File", "null"], [".bai"]),
                cwlout(["hla", "fastq"], ["File", "null"])],
               ["biobambam"],
               {"files": 3})]
    return align

def variant():
    """Variant calling workflow definition for CWL generation.
    """
    align_wf = _variant_shared()
    vc_wf = [s("get_parallel_regions", "batch-split",
               [["batch_rec"]],
               [cwlout(["region"], "string")]),
             s("variantcall_batch_region", "batch-parallel",
               [["batch_rec"], ["region"]],
               [cwlout(["vrn_file_region"], "File", [".tbi"]),
                cwlout(["region"], "string")]),
             s("concat_batch_variantcalls", "batch-merge",
               [["batch_rec"], ["region"], ["vrn_file_region"]],
               [cwlout(["vrn_file"], "File", [".tbi"])]),
             s("postprocess_variants", "batch-single",
               [["batch_rec"], ["vrn_file"]],
               [cwlout(["vrn_file"], "File", [".tbi"])]),
             s("compare_to_rm", "batch-single",
               [["batch_rec"], ["vrn_file"]],
               [cwlout(["validate", "summary"], ["File", "null"]),
                cwlout(["validate", "tp"], ["File", "null"], [".tbi"]),
                cwlout(["validate", "fp"], ["File", "null"], [".tbi"]),
                cwlout(["validate", "fn"], ["File", "null"], [".tbi"])]),
             s("vc_output_record", "batch-single",
               [["batch_rec"], ["vrn_file"], ["validate", "summary"],
                ["validate", "tp"], ["validate", "fp"], ["validate", "fn"]],
               [cwlout("vc_rec", "record")])]
    align = [w("alignment", "multi-parallel", align_wf,
               [["align_split"], ["files"], ["work_bam"], ["config", "algorithm", "quality_format"]]),
             s("prep_samples", "multi-parallel",
               [["config", "algorithm", "coverage"],
                ["config", "algorithm", "variant_regions"],
                ["reference", "fasta", "base"]],
               [cwlout(["config", "algorithm", "variant_regions"], ["File", "null"]),
                cwlout(["config", "algorithm", "variant_regions_merged"], ["File", "null"]),
                cwlout(["config", "algorithm", "variant_regions_orig"], ["File", "null"]),
                cwlout(["config", "algorithm", "coverage"], ["File", "null"]),
                cwlout(["config", "algorithm", "coverage_merged"], ["File", "null"]),
                cwlout(["config", "algorithm", "coverage_orig"], ["File", "null"]),
                cwlout(["config", "algorithm", "seq2c_bed_ready"], ["File", "null"])]),
             s("postprocess_alignment", "multi-parallel",
               [["align_bam"],
                ["config", "algorithm", "coverage_interval"],
                ["config", "algorithm", "variant_regions"],
                ["config", "algorithm", "variant_regions_merged"],
                ["config", "algorithm", "variant_regions_orig"],
                ["config", "algorithm", "coverage"],
                ["config", "algorithm", "coverage_merged"],
                ["config", "algorithm", "coverage_orig"],
                ["config", "algorithm", "seq2c_bed_ready"],
                ["config", "algorithm", "recalibrate"],
                ["reference", "fasta", "base"]],
               [cwlout(["config", "algorithm", "coverage_interval"], "string"),
                cwlout(["config", "algorithm", "variant_regions"], "File"),
                cwlout(["config", "algorithm", "variant_regions_merged"], "File"),
                cwlout(["config", "algorithm", "variant_regions_orig"], "File"),
                cwlout(["config", "algorithm", "coverage"], ["File", "null"]),
                cwlout(["config", "algorithm", "coverage_merged"], ["File", "null"]),
                cwlout(["config", "algorithm", "coverage_orig"], ["File", "null"]),
                cwlout(["config", "algorithm", "seq2c_bed_ready"], "File"),
                cwlout(["regions", "callable"], "File"),
                cwlout(["regions", "sample_callable"], "File"),
                cwlout(["regions", "nblock"], "File"),
                cwlout(["regions", "highdepth"], ["File", "null"])]),
             s("combine_sample_regions", "multi-combined",
               [["regions", "callable"], ["regions", "nblock"],
                ["config", "algorithm", "nomap_split_size"], ["config", "algorithm", "nomap_split_targets"],
                ["reference", "fasta", "base"]],
               [cwlout(["config", "algorithm", "callable_regions"], "File"),
                cwlout(["config", "algorithm", "non_callable_regions"], "File"),
                cwlout(["config", "algorithm", "callable_count"], "int")])]
    hla = [s("call_hla", "multi-parallel",
             [["hla", "fastq"]],
             [cwlout(["hla", "hlacaller"], ["string", "null"]),
              cwlout(["hla", "call_file"], ["File", "null"])])]
    vc = [s("batch_for_variantcall", "multi-batch",
            [["analysis"], ["genome_build"], ["align_bam"], ["config", "algorithm", "callable_regions"],
             ["metadata", "batch"], ["metadata", "phenotype"],
             ["regions", "callable"], ["config", "algorithm", "variantcaller"],
             ["config", "algorithm", "coverage_interval"],
             ["config", "algorithm", "variant_regions"],
             ["config", "algorithm", "validate"], ["config", "algorithm", "validate_regions"],
             ["config", "algorithm", "tools_on"],
             ["config", "algorithm", "tools_off"],
             ["reference", "fasta", "base"], ["reference", "rtg"], ["reference", "genome_context"],
             ["genome_resources", "variation", "cosmic"], ["genome_resources", "variation", "dbsnp"]],
            [cwlout("batch_rec", "record")],
            unlist=[["config", "algorithm", "variantcaller"]]),
          w("variantcall", "multi-parallel", vc_wf,
            [["region"], ["vrn_file_region"], ["vrn_file"], ["validate", "summary"]]),
          s("summarize_grading_vc", "multi-combined",
            [["vc_rec"]],
            [cwlout(["validate", "grading_summary"], ["File", "null"]),
             cwlout(["validate", "grading_plots"], {"type": "array", "items": ["File", "null"]})])]
    qc = [s("qc_to_rec", "multi-batch",
            [["align_bam"], ["analysis"], ["reference", "fasta", "base"],
             ["genome_build"], ["config", "algorithm", "coverage_interval"],
             ["config", "algorithm", "tools_on"], ["config", "algorithm", "tools_off"],
             ["config", "algorithm", "qc"],
             ["config", "algorithm", "variant_regions"],
             ["config", "algorithm", "variant_regions_merged"],
             ["config", "algorithm", "coverage"],
             ["config", "algorithm", "coverage_merged"]],
            [cwlout("qc_rec", "record")]),
          s("pipeline_summary", "multi-parallel",
            ["qc_rec"],
            [cwlout(["summary", "qc"], ["File", "null"]),
             cwlout(["summary", "metrics"], "string")]),
          s("multiqc_summary", "multi-combined",
            [["qc_rec"], ["summary", "qc"], ["summary", "metrics"]],
            [cwlout(["summary", "multiqc"], ["File", "null"])])]
    steps = align + vc + qc
    final_outputs = [["align_bam"], ["summary", "multiqc"]]
    return steps, final_outputs

def sv():
    """Structural variant workflow, for development purposes.

    Will eventually merge with variant workflow using selectors to determine
    required steps, but this is an initial step for testing and to work on
    understanding necessary steps for each process.
    """
    align = _variant_shared()
    sv = [s("detect_sv", "batch-single",
            [["sv_batch_rec"]],
            [cwlout(["sv", "0", "variantcaller"], "string"),
             cwlout(["sv", "0", "vrn_file"], "File", [".tbi"])]),
         ]
    steps = [w("alignment", "multi-parallel", align,
               [["align_split"], ["files"], ["work_bam"], ["config", "algorithm", "quality_format"]]),
             s("batch_for_sv", "multi-batch",
               [["analysis"], ["genome_build"], ["align_bam"],
                ["work_bam_plus", "disc"], ["work_bam_plus", "sr"],
                ["metadata", "batch"], ["metadata", "phenotype"],
                ["config", "algorithm", "svcaller"],
                ["config", "algorithm", "tools_on"],
                ["config", "algorithm", "tools_off"],
                ["reference", "fasta", "base"]],
               [cwlout("sv_batch_rec", "record")],
               unlist=[["config", "algorithm", "svcaller"]]),
             w("svcall", "multi-parallel", sv, []),
            ]
    final_outputs = [["align_bam"]]
    return steps, final_outputs

def fastrnaseq():
    prep = [s("prep_samples", "multi-parallel",
              [["files"],
               dd.get_keys("sample_name")],
              [cwlout(["files"], "File")],
              programs=["picard"])]
    quant = [s("run_salmon_reads", "multi-parallel",
               [["files"],
                dd.get_keys("sample_name"),
                dd.get_keys("gtf_file"),
                dd.get_keys("ref_file"),
                dd.get_keys("genome_build")],
               [cwlout(dd.get_keys("sailfish_dir"), "File")],
               programs=["salmon"],
               disk={"files": 1.5})]
    steps = quant
    final_outputs = [dd.get_keys('sailfish_dir')]
    return steps, final_outputs

def rnaseq():
    prep = [s("prep_samples", "multi-parallel",
              [["files"],
               dd.get_keys("sample_name")],
              [cwlout(["files"], "File")],
              programs=["picard"])]
    align = [s("process_alignment", "multi-parallel",
               [["files"], ["reference", "fasta", "base"],
                ["analysis"],
                ["rgnames", "pl"], ["rgnames", "sample"], ["rgnames", "pu"],
                ["rgnames", "lane"], ["rgnames", "rg"], ["rgnames", "lb"],
                ["reference", "aligner", "indexes"],
                ["config", "algorithm", "aligner"],
                ["genome_resources", "rnaseq", "transcripts"],
                ["config", "algorithm", "quality_format"]],
               [cwlout(["work_bam"], "File"),
                cwlout(["align_bam"], "File")],
               ["aligner", "samtools", "sambamba"],
               {"files": 1.5})]
    quantitate = [s("rnaseq_quantitate", "multi-parallel",
                  [["files"],
                   dd.get_keys("work_bam"),
                   dd.get_keys("gtf_file"),
                   dd.get_keys("ref_file"),
                   dd.get_keys("genome_build")],
                  [cwlout(dd.get_keys("count_file"), "File"),
                   cwlout(dd.get_keys("sailfish_dir"), "File")],
                  programs=["sailfish"],
                  disk={"files": 1.5})]
    qc = [s("pipeline_summary", "multi-parallel",
            [["align_bam"], ["analysis"], ["reference", "fasta", "base"],
             ["config", "algorithm", "qc"]],
            [cwlout(["summary", "qc", "samtools"], "File"),
             cwlout(["summary", "qc", "fastqc"], "File")],
            ["samtools", "fastqc"]),
          s("multiqc_summary", "multi-combined",
            [["genome_build"], ["summary", "qc", "samtools"], ["summary", "qc", "fastqc"],
             ["reference", "fasta", "base"], ["config", "algorithm", "coverage_interval"]],
            [cwlout(["summary", "multiqc"], ["File", "null"])])]

    steps = prep + align + quantitate + qc
    final_outputs = [dd.get_keys("work_bam"), dd.get_keys("sailfish_dir"),
                     ["summary", "multiqc"]]
    return steps, final_outputs

workflows = \
  {"variant": variant, "variant2": variant, "fastrna-seq": fastrnaseq,
   "rna-seq": rnaseq, "sv": sv}
