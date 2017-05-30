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

def s(name, parallel, inputs, outputs, image, programs=None, disk=None, cores=None, unlist=None):
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
    Step = collections.namedtuple("Step", "name parallel inputs outputs image programs disk cores unlist")
    if programs is None: programs = []
    if unlist is None: unlist = []
    return Step(name, parallel, inputs, outputs, image, programs, disk, cores, unlist)

def w(name, parallel, workflow, internal):
    """A workflow, allowing specification of sub-workflows for nested parallelization.

    name and parallel are documented under the Step (s) function.
    workflow -- a list of Step tuples defining the sub-workflow
    internal -- variables used in the sub-workflow but not exposed to subsequent steps
    """
    Workflow = collections.namedtuple("Workflow", "name parallel workflow internal")
    return Workflow(name, parallel, workflow, internal)

def cwlout(key, valtype=None, extensions=None, fields=None):
    """Definition of an output variable, defining the type and associated secondary files.
    """
    out = {"id": key}
    if valtype:
        out["type"] = valtype
    if fields:
        out["fields"] = fields
    if extensions:
        out["secondaryFiles"] = extensions
    return out

def _variant_shared():
    align = [s("prep_align_inputs", "single-split",
               [["alignment_rec"]],
               [cwlout("process_alignment_rec", "record",
                       fields=[cwlout(["files"], {"type": "array", "items": "File"}),
                               cwlout(["config", "algorithm", "quality_format"], "string"),
                               cwlout(["align_split"], ["string", "null"])])],
               "bcbio-align", ["grabix", "htslib", "biobambam"],
               {"files": 1.5}, cores=1),
             s("process_alignment", "single-parallel",
               [["alignment_rec"], ["process_alignment_rec"]],
               [cwlout(["work_bam"], "File"),
                cwlout(["align_bam"], "File"),
                cwlout(["hla", "fastq"], ["File", "null", {"type": "array", "items": "File"}]),
                cwlout(["work_bam_plus", "disc"], ["File", "null"], [".bai"]),
                cwlout(["work_bam_plus", "sr"], ["File", "null"], [".bai"])],
               "bcbio-align", ["bwa", "bwakit", "grabix", "novoalign", "snap-aligner=1.0dev.97",
                               "sentieon", "samtools", "sambamba", "fgbio", "umis", "biobambam", "seqtk",
                               "samblaster"],
               {"files": 1.5}),
             s("merge_split_alignments", "single-merge",
               [["alignment_rec"], ["work_bam"], ["align_bam"],
                ["work_bam_plus", "disc"], ["work_bam_plus", "sr"],["hla", "fastq"]],
               [cwlout(["align_bam"], "File", [".bai"]),
                cwlout(["work_bam_plus", "disc"], ["File", "null"], [".bai"]),
                cwlout(["work_bam_plus", "sr"], ["File", "null"], [".bai"]),
                cwlout(["hla", "fastq"], ["File", "null", {"type": "array", "items": "File"}])],
               "bcbio-align", ["biobambam", "samtools"],
               {"files": 3})]
    return align

def _variant_hla(checkpoints):
    """Add hla analysis to workflow, if configured.
    """
    if not checkpoints.get("hla"):
        return [], []
    hla = [s("call_hla", "multi-parallel",
             [["hla", "fastq"]],
             [cwlout(["hla", "hlacaller"], ["string", "null"]),
              cwlout(["hla", "call_file"], ["File", "null"])],
             "bcbio-vc", ["optitype"])]
    return hla, []

def _variant_vc(checkpoints):
    """Add variant calling to workflow, if configured.
    """
    if not checkpoints.get("vc"):
        return [], []
    vc_wf = [s("get_parallel_regions", "batch-split",
               [["batch_rec"]],
               [cwlout(["region"], "string")],
               "bcbio-base",
               cores=1),
             s("variantcall_batch_region", "batch-parallel",
               [["batch_rec"], ["region"]],
               [cwlout(["vrn_file_region"], "File", [".tbi"]),
                cwlout(["region"], "string")],
               "bcbio-vc", ["bcftools", "bedtools", "freebayes=1.1.0",
                            "gatk", "gatk4", "gatk-framework",
                            "htslib", "picard", "platypus-variant", "pythonpy",
                            "samtools", "vardict", "vardict-java", "varscan", "vcfanno",
                            "vcflib", "vt", "r=3.2.2", "perl"],
               cores=1),
             s("concat_batch_variantcalls", "batch-merge",
               [["batch_rec"], ["region"], ["vrn_file_region"]],
               [cwlout(["vrn_file"], "File", [".tbi"])],
               "bcbio-vc", ["bcftools", "htslib"],
               cores=1),
             s("postprocess_variants", "batch-single",
               [["batch_rec"], ["vrn_file"]],
               [cwlout(["vrn_file"], "File", [".tbi"])],
               "bcbio-vc", ["snpeff=4.3i"]),
             s("compare_to_rm", "batch-single",
               [["batch_rec"], ["vrn_file"]],
               [cwlout("vc_rec", "record",
                       fields=[cwlout(["validate", "summary"], ["File", "null"]),
                               cwlout(["validate", "tp"], ["File", "null"], [".tbi"]),
                               cwlout(["validate", "fp"], ["File", "null"], [".tbi"]),
                               cwlout(["validate", "fn"], ["File", "null"], [".tbi"]),
                               cwlout("inherit")])],
               "bcbio-vc", ["bcftools", "bedtools", "pythonpy", "gvcf-regions",
                            "htslib", "rtg-tools", "vcfanno"])]
    vc = [s("batch_for_variantcall", "multi-batch",
            [["analysis"], ["genome_build"], ["align_bam"], ["config", "algorithm", "callable_regions"],
             ["metadata", "batch"], ["metadata", "phenotype"],
             ["regions", "sample_callable"], ["config", "algorithm", "variantcaller"],
             ["config", "algorithm", "coverage_interval"],
             ["config", "algorithm", "variant_regions"],
             ["config", "algorithm", "validate"], ["config", "algorithm", "validate_regions"],
             ["config", "algorithm", "tools_on"],
             ["config", "algorithm", "tools_off"],
             ["reference", "fasta", "base"], ["reference", "rtg"], ["reference", "genome_context"],
             ["genome_resources", "variation", "cosmic"], ["genome_resources", "variation", "dbsnp"]],
            [cwlout("batch_rec", "record")],
            "bcbio-base",
            cores=1,
            unlist=[["config", "algorithm", "variantcaller"]]),
          w("variantcall", "multi-parallel", vc_wf,
            [["region"], ["vrn_file_region"], ["vrn_file"], ["validate", "summary"]]),
          s("summarize_grading_vc", "multi-combined",
            [["vc_rec"]],
            [cwlout(["validate", "grading_summary"], ["File", "null"]),
             cwlout(["validate", "grading_plots"], {"type": "array", "items": ["File", "null"]})],
            "bcbio-base",
            cores=1)]
    return vc, [["validate", "grading_summary"]]

def _variant_checkpoints(samples):
    """Check sample configuration to identify required steps in analysis.
    """
    checkpoints = {}
    checkpoints["vc"] = any([dd.get_variantcaller(d) for d in samples])
    checkpoints["hla"] = any([dd.get_hlacaller(d) for d in samples])
    return checkpoints

def variant(samples):
    """Variant calling workflow definition for CWL generation.
    """
    checkpoints = _variant_checkpoints(samples)
    align_wf = _variant_shared()
    align = [s("alignment_to_rec", "multi-combined",
               [["files"],
                ["config", "algorithm", "align_split_size"],
                ["reference", "fasta", "base"],
                ["rgnames", "pl"], ["rgnames", "sample"], ["rgnames", "pu"],
                ["rgnames", "lane"], ["rgnames", "rg"], ["rgnames", "lb"],
                ["reference", "aligner", "indexes"],
                ["config", "algorithm", "aligner"],
                ["config", "algorithm", "mark_duplicates"]],
               [cwlout("alignment_rec", "record")],
               "bcbio-base",
               cores=1),
             w("alignment", "multi-parallel", align_wf,
               [["align_split"], ["process_alignment_rec"],
                ["work_bam"], ["config", "algorithm", "quality_format"]]),
             s("prep_samples_to_rec", "multi-combined",
               [["config", "algorithm", "coverage"],
                ["config", "algorithm", "variant_regions"],
                ["reference", "fasta", "base"]],
               [cwlout("prep_samples_rec", "record")],
               "bcbio-base",
               cores=1),
             s("prep_samples", "multi-parallel",
               ["prep_samples_rec"],
               [cwlout(["config", "algorithm", "variant_regions"], ["File", "null"]),
                cwlout(["config", "algorithm", "variant_regions_merged"], ["File", "null"]),
                cwlout(["config", "algorithm", "variant_regions_orig"], ["File", "null"]),
                cwlout(["config", "algorithm", "coverage"], ["File", "null"]),
                cwlout(["config", "algorithm", "coverage_merged"], ["File", "null"]),
                cwlout(["config", "algorithm", "coverage_orig"], ["File", "null"]),
                cwlout(["config", "algorithm", "seq2c_bed_ready"], ["File", "null"])],
               "bcbio-align", ["htslib", "bedtools", "pythonpy"],
               cores=1),
             s("postprocess_alignment_to_rec", "multi-combined",
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
               [cwlout("postprocess_alignment_rec", "record")],
               "bcbio-base",
               cores=1),
             s("postprocess_alignment", "multi-parallel",
               [["postprocess_alignment_rec"]],
               [cwlout(["config", "algorithm", "coverage_interval"], "string"),
                cwlout(["config", "algorithm", "variant_regions"], "File"),
                cwlout(["config", "algorithm", "variant_regions_merged"], "File"),
                cwlout(["config", "algorithm", "variant_regions_orig"], "File"),
                cwlout(["config", "algorithm", "coverage"], ["File", "null"]),
                cwlout(["config", "algorithm", "coverage_merged"], ["File", "null"]),
                cwlout(["config", "algorithm", "coverage_orig"], ["File", "null"]),
                cwlout(["config", "algorithm", "seq2c_bed_ready"], ["File", "null"]),
                cwlout(["regions", "callable"], "File"),
                cwlout(["regions", "sample_callable"], "File"),
                cwlout(["regions", "nblock"], "File"),
                cwlout(["regions", "highdepth"], ["File", "null"])],
               "bcbio-align", ["sambamba", "goleft", "bedtools", "htslib"]),
             s("combine_sample_regions", "multi-combined",
               [["regions", "callable"], ["regions", "nblock"],
                ["config", "algorithm", "nomap_split_size"], ["config", "algorithm", "nomap_split_targets"],
                ["reference", "fasta", "base"]],
               [cwlout(["config", "algorithm", "callable_regions"], "File"),
                cwlout(["config", "algorithm", "non_callable_regions"], "File"),
                cwlout(["config", "algorithm", "callable_count"], "int")],
               "bcbio-align", ["bedtools", "htslib"],
               cores=1)]
    qc = [s("qc_to_rec", "multi-combined",
            [["align_bam"], ["analysis"], ["reference", "fasta", "base"],
             ["genome_build"], ["config", "algorithm", "coverage_interval"],
             ["config", "algorithm", "tools_on"], ["config", "algorithm", "tools_off"],
             ["config", "algorithm", "qc"],
             ["config", "algorithm", "variant_regions"],
             ["config", "algorithm", "variant_regions_merged"],
             ["config", "algorithm", "coverage"],
             ["config", "algorithm", "coverage_merged"]],
            [cwlout("qc_rec", "record")],
            "bcbio-base",
            cores=1),
          s("pipeline_summary", "multi-parallel",
            ["qc_rec"],
            [cwlout("qcout_rec", "record",
                    fields=[cwlout(["summary", "qc"], ["File", "null"]),
                            cwlout(["summary", "metrics"], "string"),
                            cwlout("inherit")])],
            "bcbio-qc", ["bcftools", "bedtools", "fastqc", "goleft", "picard", "pythonpy",
                         "qsignature", "qualimap", "sambamba", "samtools"]),
          s("multiqc_summary", "multi-combined",
            [["qcout_rec"]],
            [cwlout(["summary", "multiqc"], ["File", "null"])],
            "bcbio-qc", ["multiqc", "multiqc-bcbio"],
            cores=1)]
    vc, vc_out = _variant_vc(checkpoints)
    hla, hla_out = _variant_hla(checkpoints)
    steps = align + qc + vc + hla
    final_outputs = [["align_bam"],
                     cwlout(["summary", "multiqc"], {"type": "array", "items": ["File", "null"]}),
                    ] + vc_out + hla_out
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
             cwlout(["sv", "0", "vrn_file"], "File", [".tbi"])],
            "bcbio"),
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
               "bcbio-base",
               unlist=[["config", "algorithm", "svcaller"]]),
             w("svcall", "multi-parallel", sv, []),
            ]
    final_outputs = [["align_bam"]]
    return steps, final_outputs

def fastrnaseq(samples):
    prep = [s("prep_samples", "multi-parallel",
              [["files"],
               dd.get_keys("sample_name")],
              [cwlout(["files"], "File")],
              "bcbio", programs=["picard"])]
    quant = [s("run_salmon_reads", "multi-parallel",
               [["files"],
                dd.get_keys("sample_name"),
                dd.get_keys("gtf_file"),
                dd.get_keys("ref_file"),
                dd.get_keys("genome_build")],
               [cwlout(dd.get_keys("sailfish_dir"), "File")],
               "bcbio", programs=["salmon"],
               disk={"files": 1.5})]
    steps = quant
    final_outputs = [dd.get_keys('sailfish_dir')]
    return steps, final_outputs

def rnaseq(samples):
    prep = [s("prep_samples", "multi-parallel",
              [["files"],
               dd.get_keys("sample_name")],
              [cwlout(["files"], "File")],
              "bcbio", programs=["picard"])]
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
               "bcbio-align", ["aligner", "samtools", "sambamba", "seqtk"],
               {"files": 1.5})]
    quantitate = [s("rnaseq_quantitate", "multi-parallel",
                  [["files"],
                   dd.get_keys("work_bam"),
                   dd.get_keys("gtf_file"),
                   dd.get_keys("ref_file"),
                   dd.get_keys("genome_build")],
                  [cwlout(dd.get_keys("count_file"), "File"),
                   cwlout(dd.get_keys("sailfish_dir"), "File")],
                  "bcbio", programs=["sailfish"],
                  disk={"files": 1.5})]
    qc = [s("pipeline_summary", "multi-parallel",
            [["align_bam"], ["analysis"], ["reference", "fasta", "base"],
             ["config", "algorithm", "qc"]],
            [cwlout(["summary", "qc", "samtools"], "File"),
             cwlout(["summary", "qc", "fastqc"], "File")],
            "bcbio", ["samtools", "fastqc"]),
          s("multiqc_summary", "multi-combined",
            [["genome_build"], ["summary", "qc", "samtools"], ["summary", "qc", "fastqc"],
             ["reference", "fasta", "base"], ["config", "algorithm", "coverage_interval"]],
            [cwlout(["summary", "multiqc"], ["File", "null"])],
            "bcbio")]

    steps = prep + align + quantitate + qc
    final_outputs = [dd.get_keys("work_bam"), dd.get_keys("sailfish_dir"),
                     ["summary", "multiqc"]]
    return steps, final_outputs

workflows = \
  {"variant": variant, "variant2": variant, "fastrna-seq": fastrnaseq,
   "rna-seq": rnaseq, "sv": sv}
