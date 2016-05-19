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

def s(name, parallel, inputs, outputs, programs=None, disk=None, noinputs=None, unlist=None):
    """Represent a step in a workflow.

    name -- The run function name, which must match a definition in distributed/multitasks
    inputs -- List of input keys required for the function. Each key is of the type:
      ["toplevel", "sublevel"] -- an argument you could pass to toolz.get_in.
      You only define file based inputs, non-file inputs are always passed unless
      specified in noinputs.
    outputs -- List of outputs with information about file type. Use cwlout functions
    programs -- Required programs for this step, used to define resource usage.
    disk -- Information about disk usage requirements, specified as multipliers of
            input files. Ensures enough disk present when that is a limiting factor
            when selecting cloud node resources.
    noinputs -- Non-file variables to not pass to this step
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
    Step = collections.namedtuple("Step", "name parallel inputs outputs programs disk noinputs unlist")
    if programs is None: programs = []
    if noinputs is None: noinputs = []
    if unlist is None: unlist = []
    return Step(name, parallel, inputs, outputs, programs, disk, noinputs, unlist)

def w(name, parallel, workflow, internal, noinputs=None):
    """A workflow, allowing specification of sub-workflows for nested parallelization.

    name, parallel and noinputs are documented under the Step (s) function.
    workflow -- a list of Step tuples defining the sub-workflow
    internal -- variables used in the sub-workflow but not exposed to subsequent steps
    """
    Workflow = collections.namedtuple("Workflow", "name parallel workflow internal noinputs")
    if noinputs is None: noinputs = []
    return Workflow(name, parallel, workflow, internal, noinputs)

def cwlout(key, valtype, extensions=None):
    """Definition of an output variable, defining the type and associated secondary files.
    """
    out = {"id": key,
           "type": valtype}
    if extensions:
        out["secondaryFiles"] = extensions
    return out

def variant():
    """Variant calling workflow definition for CWL generation.
    """
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
               [["files"], ["reference", "fasta", "base"],
                ["rgnames", "pl"], ["rgnames", "sample"], ["rgnames", "pu"],
                ["rgnames", "lane"], ["rgnames", "rg"], ["rgnames", "lb"],
                ["reference", "aligner", "indexes"],
                ["config", "algorithm", "aligner"],
                ["config", "algorithm", "mark_duplicates"]],
               [cwlout(["work_bam"], "File"),
                cwlout(["align_bam"], "File"),
                cwlout(["hla", "fastq"], ["File", "null"]),
                cwlout(["work_bam-plus", "disc"], "File", [".bai"]),
                cwlout(["work_bam-plus", "sr"], "File", [".bai"])],
               ["aligner", "samtools", "sambamba"],
               {"files": 1.5}),
             s("merge_split_alignments", "single-merge",
               [["work_bam"], ["align_bam"], ["work_bam-plus", "disc"], ["work_bam-plus", "sr"],
                ["hla", "fastq"]],
               [cwlout(["align_bam"], "File", [".bai"]),
                cwlout(["work_bam-plus", "disc"], "File", [".bai"]),
                cwlout(["work_bam-plus", "sr"], "File", [".bai"]),
                cwlout(["hla", "fastq"], ["File", "null"])],
               ["biobambam"],
               {"files": 3},
               noinputs=[["align_split"], ["config", "algorithm", "quality_format"]])]
    vc = [s("get_parallel_regions", "batch-split",
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
            [cwlout(["vrn_file"], "File", [".tbi"])],
            noinputs=[["region"]]),
          s("compare_to_rm", "batch-single",
            [["batch_rec"], ["vrn_file"]],
            [cwlout(["validate", "summary"], ["File", "null"]),
             cwlout(["validate", "tp"], ["File", "null"], [".tbi"]),
             cwlout(["validate", "fp"], ["File", "null"], [".tbi"]),
             cwlout(["validate", "fn"], ["File", "null"], [".tbi"])],
            noinputs=[["region"]])]
    steps = [w("alignment", "multi-parallel", align,
               [["align_split"], ["files"], ["work_bam"], ["config", "algorithm", "quality_format"]]),
             s("prep_samples", "multi-parallel",
               [["config", "algorithm", "variant_regions"],
                ["reference", "fasta", "base"]],
               [cwlout(["config", "algorithm", "variant_regions"], ["File", "null"]),
                cwlout(["config", "algorithm", "variant_regions_merged"], ["File", "null"])]),
             s("postprocess_alignment", "multi-parallel",
               [["align_bam"], ["config", "algorithm", "variant_regions"],
                ["config", "algorithm", "variant_regions_merged"],
                ["config", "algorithm", "recalibrate"],
                ["reference", "fasta", "base"]],
               [cwlout(["config", "algorithm", "coverage_interval"], "string"),
                cwlout(["config", "algorithm", "variant_regions"], "File"),
                cwlout(["config", "algorithm", "variant_regions_merged"], "File"),
                cwlout(["regions", "callable"], "File"),
                cwlout(["regions", "sample_callable"], "File"),
                cwlout(["regions", "nblock"], "File"),
                cwlout(["regions", "highdepth"], ["File", "null"]),
                cwlout(["regions", "offtarget_stats"], ["File", "null"])]),
             # s("call_hla", "multi-parallel",
             #   [["hla", "fastq"]],
             #   [cwlout(["hla", "hlacaller"], ["string", "null"]),
             #    cwlout(["hla", "call_file"], ["File", "null"])]),
             s("combine_sample_regions", "multi-combined",
               [["regions", "callable"], ["regions", "nblock"],
                ["config", "algorithm", "nomap_split_size"], ["config", "algorithm", "nomap_split_targets"],
                ["reference", "fasta", "base"]],
               [cwlout(["config", "algorithm", "callable_regions"], "File"),
                cwlout(["config", "algorithm", "non_callable_regions"], "File"),
                cwlout(["config", "algorithm", "callable_count"], "int")]),
             s("batch_for_variantcall", "multi-batch",
               [["analysis"], ["genome_build"], ["align_bam"], ["config", "algorithm", "callable_regions"],
                #["metadata", "batch"], ["metadata", "phenotype"],
                ["regions", "callable"], ["config", "algorithm", "variantcaller"],
                ["config", "algorithm", "variant_regions"],
                ["config", "algorithm", "validate"], ["config", "algorithm", "validate_regions"],
                ["config", "algorithm", "tools_on"],
                ["reference", "fasta", "base"], ["reference", "rtg"],
                ["genome_resources", "variation", "cosmic"], ["genome_resources", "variation", "dbsnp"]],
               "batch_rec",
               unlist=[["config", "algorithm", "variantcaller"]],
               noinputs=[["hla", "hlacaller"], ["config", "algorithm", "callable_count"]]),
             w("variantcall", "multi-parallel", vc,
               [["region"], ["vrn_file_region"]],
               noinputs=[["hla", "hlacaller"], ["config", "algorithm", "callable_count"]]),
             s("pipeline_summary", "multi-parallel",
               [["align_bam"], ["analysis"], ["reference", "fasta", "base"],
                ["config", "algorithm", "coverage"],
                ["config", "algorithm", "qc"],
                ["config", "algorithm", "variant_regions"], ["regions", "offtarget_stats"]],
               [cwlout(["summary", "qc", "samtools"], "File"),
                cwlout(["summary", "qc", "fastqc"], "File"),
                cwlout(["coverage", "all"], ["File", "null"]),
                cwlout(["coverage", "problems"], ["File", "null"])],
               ["samtools", "fastqc"]),
             s("multiqc_summary", "multi-combined",
               [["summary", "qc", "samtools"], ["summary", "qc", "fastqc"]],
               [cwlout(["summary", "multiqc"], ["File", "null"])])
             # s("qc_report_summary", "multi-combined",
             #   [["align_bam"],
             #    ["reference", "fasta", "base"],
             #    ["summary", "qc"], ["coverage", "all"], ["coverage", "problems"]],
             #   [cwlout(["coverage", "report"], ["File", "null"])])
             ]
    final_outputs = [["align_bam"], ["vrn_file"], ["validate", "summary"], ["summary", "multiqc"]]
    return steps, final_outputs

workflows = \
  {"variant": variant, "variant2": variant}
