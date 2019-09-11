"""Main entry point for distributed next-gen sequencing pipelines.

Handles running the full pipeline based on instructions
"""
from __future__ import print_function
from collections import defaultdict
import copy
import os
import sys
import resource
import tempfile

import toolz as tz

from bcbio import log, heterogeneity, hla, structural, utils
from bcbio.cwl.inspect import initialize_watcher
from bcbio.distributed import prun
from bcbio.distributed.transaction import tx_tmpdir
from bcbio.log import logger, DEFAULT_LOG_DIR
from bcbio.ngsalign import alignprep
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import (archive, config_utils, disambiguate, region,
                            run_info, qcsummary, rnaseq)
from bcbio.provenance import profile, system
from bcbio.variation import (ensemble, genotype, population, validate, joint,
                             peddy)
from bcbio.chipseq import peaks

def run_main(workdir, config_file=None, fc_dir=None, run_info_yaml=None,
             parallel=None, workflow=None):
    """Run variant analysis, handling command line options.
    """
    # Set environment to standard to use periods for decimals and avoid localization
    locale_to_use = utils.get_locale()
    os.environ["LC_ALL"] = locale_to_use
    os.environ["LC"] = locale_to_use
    os.environ["LANG"] = locale_to_use
    workdir = utils.safe_makedir(os.path.abspath(workdir))
    os.chdir(workdir)
    config, config_file = config_utils.load_system_config(config_file, workdir)
    parallel = log.create_base_logger(config, parallel)
    log.setup_local_logging(config, parallel)
    logger.info(f"System YAML configuration: {os.path.abspath(config_file)}.")
    logger.info(f"Locale set to {locale_to_use}.")
    if config.get("log_dir", None) is None:
        config["log_dir"] = os.path.join(workdir, DEFAULT_LOG_DIR)
    if parallel["type"] in ["local", "clusterk"]:
        _setup_resources()
        _run_toplevel(config, config_file, workdir, parallel,
                      fc_dir, run_info_yaml)
    elif parallel["type"] == "ipython":
        assert parallel["scheduler"] is not None, "IPython parallel requires a specified scheduler (-s)"
        if parallel["scheduler"] != "sge":
            assert parallel["queue"] is not None, "IPython parallel requires a specified queue (-q)"
        elif not parallel["queue"]:
            parallel["queue"] = ""
        _run_toplevel(config, config_file, workdir, parallel,
                      fc_dir, run_info_yaml)
    else:
        raise ValueError("Unexpected type of parallel run: %s" % parallel["type"])

def _setup_resources():
    """Attempt to increase resource limits up to hard limits.

    This allows us to avoid out of file handle limits where we can
    move beyond the soft limit up to the hard limit.
    """
    target_procs = 10240
    cur_proc, max_proc = resource.getrlimit(resource.RLIMIT_NPROC)
    target_proc = min(max_proc, target_procs) if max_proc > 0 else target_procs
    resource.setrlimit(resource.RLIMIT_NPROC, (max(cur_proc, target_proc), max_proc))
    cur_hdls, max_hdls = resource.getrlimit(resource.RLIMIT_NOFILE)
    target_hdls = min(max_hdls, target_procs) if max_hdls > 0 else target_procs
    resource.setrlimit(resource.RLIMIT_NOFILE, (max(cur_hdls, target_hdls), max_hdls))

def _run_toplevel(config, config_file, work_dir, parallel,
                  fc_dir=None, run_info_yaml=None):
    """
    Run toplevel analysis, processing a set of input files.
    config_file -- Main YAML configuration file with system parameters
    fc_dir -- Directory of fastq files to process
    run_info_yaml -- YAML configuration file specifying inputs to process
    """
    dirs = run_info.setup_directories(work_dir, fc_dir, config, config_file)
    config_file = os.path.join(dirs["config"], os.path.basename(config_file))
    pipelines, config = _pair_samples_with_pipelines(run_info_yaml, config)
    system.write_info(dirs, parallel, config)
    with tx_tmpdir(config if parallel.get("type") == "local" else None) as tmpdir:
        tempfile.tempdir = tmpdir
        for pipeline, samples in pipelines.items():
            for xs in pipeline(config, run_info_yaml, parallel, dirs, samples):
                pass

# ## Generic pipeline framework

def _wres(parallel, progs, fresources=None, ensure_mem=None):
    """Add resource information to the parallel environment on required programs and files.

    Enables spinning up required machines and operating in non-shared filesystem
    environments.

    progs -- Third party tools used in processing
    fresources -- Required file-based resources needed. These will be transferred on non-shared
                  filesystems.
    ensure_mem -- Dictionary of required minimum memory for programs used. Ensures
                  enough memory gets allocated on low-core machines.
    """
    parallel = copy.deepcopy(parallel)
    parallel["progs"] = progs
    if fresources:
        parallel["fresources"] = fresources
    if ensure_mem:
        parallel["ensure_mem"] = ensure_mem
    return parallel

def variant2pipeline(config, run_info_yaml, parallel, dirs, samples):
    ## Alignment and preparation requiring the entire input file (multicore cluster)
    # Assign GATK supplied memory if required for post-process recalibration
    align_programs = ["aligner", "samtools", "sambamba"]
    if any(tz.get_in(["algorithm", "recalibrate"], utils.to_single_data(d)) in [True, "gatk"] for d in samples):
        align_programs.append("gatk")
    with prun.start(_wres(parallel, align_programs,
                            (["reference", "fasta"], ["reference", "aligner"], ["files"])),
                    samples, config, dirs, "multicore",
                    multiplier=alignprep.parallel_multiplier(samples)) as run_parallel:
        with profile.report("organize samples", dirs):
            samples = run_parallel("organize_samples", [[dirs, config, run_info_yaml,
                                                            [x[0]["description"] for x in samples]]])
        with profile.report("alignment preparation", dirs):
            samples = run_parallel("prep_align_inputs", samples)
            samples = run_parallel("disambiguate_split", [samples])
        with profile.report("alignment", dirs):
            samples = run_parallel("process_alignment", samples)
            samples = disambiguate.resolve(samples, run_parallel)
            samples = alignprep.merge_split_alignments(samples, run_parallel)
        with profile.report("callable regions", dirs):
            samples = run_parallel("prep_samples", [samples])
            samples = run_parallel("postprocess_alignment", samples)
            samples = run_parallel("combine_sample_regions", [samples])
            samples = run_parallel("calculate_sv_bins", [samples])
            samples = run_parallel("calculate_sv_coverage", samples)
            samples = run_parallel("normalize_sv_coverage", [samples])
            samples = region.clean_sample_data(samples)
        with profile.report("hla typing", dirs):
            samples = hla.run(samples, run_parallel)

    ## Variant calling on sub-regions of the input file (full cluster)
    with prun.start(_wres(parallel, ["gatk", "picard", "variantcaller"]),
                    samples, config, dirs, "full",
                    multiplier=region.get_max_counts(samples), max_multicore=1) as run_parallel:
        with profile.report("alignment post-processing", dirs):
            samples = region.parallel_prep_region(samples, run_parallel)
        with profile.report("variant calling", dirs):
            samples = genotype.parallel_variantcall_region(samples, run_parallel)

    ## Finalize variants, BAMs and population databases (per-sample multicore cluster)
    with prun.start(_wres(parallel, ["gatk", "gatk-vqsr", "snpeff", "bcbio_variation",
                                     "gemini", "samtools", "fastqc", "sambamba",
                                     "bcbio-variation-recall", "qsignature",
                                     "svcaller", "kraken", "preseq"]),
                    samples, config, dirs, "multicore2",
                    multiplier=structural.parallel_multiplier(samples)) as run_parallel:
        with profile.report("joint squaring off/backfilling", dirs):
            samples = joint.square_off(samples, run_parallel)
        with profile.report("variant post-processing", dirs):
            samples = run_parallel("postprocess_variants", samples)
            samples = run_parallel("split_variants_by_sample", samples)
        with profile.report("prepped BAM merging", dirs):
            samples = region.delayed_bamprep_merge(samples, run_parallel)
        with profile.report("validation", dirs):
            samples = run_parallel("compare_to_rm", samples)
            samples = genotype.combine_multiple_callers(samples)
        with profile.report("ensemble calling", dirs):
            samples = ensemble.combine_calls_parallel(samples, run_parallel)
        with profile.report("validation summary", dirs):
            samples = validate.summarize_grading(samples)
        with profile.report("structural variation", dirs):
            samples = structural.run(samples, run_parallel, "initial")
        with profile.report("structural variation", dirs):
            samples = structural.run(samples, run_parallel, "standard")
        with profile.report("structural variation ensemble", dirs):
            samples = structural.run(samples, run_parallel, "ensemble")
        with profile.report("structural variation validation", dirs):
            samples = run_parallel("validate_sv", samples)
        with profile.report("heterogeneity", dirs):
            samples = heterogeneity.run(samples, run_parallel)
        with profile.report("population database", dirs):
            samples = population.prep_db_parallel(samples, run_parallel)
        with profile.report("peddy check", dirs):
            samples = peddy.run_peddy_parallel(samples, run_parallel)
        with profile.report("quality control", dirs):
            samples = qcsummary.generate_parallel(samples, run_parallel)
        with profile.report("archive", dirs):
            samples = archive.compress(samples, run_parallel)
        with profile.report("upload", dirs):
            samples = run_parallel("upload_samples", samples)
            for sample in samples:
                run_parallel("upload_samples_project", [sample])
    logger.info("Timing: finished")
    return samples

def _debug_samples(i, samples):
    print("---", i, len(samples))
    for sample in (utils.to_single_data(x) for x in samples):
        print("  ", sample["description"], sample.get("region"), \
            utils.get_in(sample, ("config", "algorithm", "variantcaller")), \
            utils.get_in(sample, ("config", "algorithm", "jointcaller")), \
            utils.get_in(sample, ("metadata", "batch")), \
            [x.get("variantcaller") for x in sample.get("variants", [])], \
            sample.get("work_bam"), \
            sample.get("vrn_file"))

def standardpipeline(config, run_info_yaml, parallel, dirs, samples):
    ## Alignment and preparation requiring the entire input file (multicore cluster)
    with prun.start(_wres(parallel, ["aligner", "samtools", "sambamba"]),
                    samples, config, dirs, "multicore") as run_parallel:
        with profile.report("organize samples", dirs):
            samples = run_parallel("organize_samples", [[dirs, config, run_info_yaml,
                                                            [x[0]["description"] for x in samples]]])
        with profile.report("alignment", dirs):
            samples = run_parallel("process_alignment", samples)
        with profile.report("callable regions", dirs):
            samples = run_parallel("prep_samples", [samples])
            samples = run_parallel("postprocess_alignment", samples)
            samples = run_parallel("combine_sample_regions", [samples])
            samples = region.clean_sample_data(samples)
    ## Quality control
    with prun.start(_wres(parallel, ["fastqc", "qsignature", "kraken", "gatk", "samtools", "preseq"]),
                    samples, config, dirs, "multicore2") as run_parallel:
        with profile.report("quality control", dirs):
            samples = qcsummary.generate_parallel(samples, run_parallel)
        with profile.report("upload", dirs):
            samples = run_parallel("upload_samples", samples)
            for sample in samples:
                run_parallel("upload_samples_project", [sample])
    logger.info("Timing: finished")
    return samples

def rnaseqpipeline(config, run_info_yaml, parallel, dirs, samples):
    samples = rnaseq_prep_samples(config, run_info_yaml, parallel, dirs, samples)
    with prun.start(_wres(parallel, ["aligner", "picard", "samtools"],
                            ensure_mem={"tophat": 10, "tophat2": 10, "star": 2, "hisat2": 8}),
                    samples, config, dirs, "alignment",
                    multiplier=alignprep.parallel_multiplier(samples)) as run_parallel:
        with profile.report("alignment", dirs):
            samples = run_parallel("disambiguate_split", [samples])
            samples = run_parallel("process_alignment", samples)
    with prun.start(_wres(parallel, ["samtools", "cufflinks"]),
                    samples, config, dirs, "rnaseqcount") as run_parallel:
        with profile.report("disambiguation", dirs):
            samples = disambiguate.resolve(samples, run_parallel)
        with profile.report("transcript assembly", dirs):
            samples = rnaseq.assemble_transcripts(run_parallel, samples)
        with profile.report("estimate expression (threaded)", dirs):
            samples = rnaseq.quantitate_expression_parallel(samples, run_parallel)

    with prun.start(_wres(parallel, ["dexseq", "express"]), samples, config,
                    dirs, "rnaseqcount-singlethread", max_multicore=1) as run_parallel:
        with profile.report("estimate expression (single threaded)", dirs):
            samples = rnaseq.quantitate_expression_noparallel(samples, run_parallel)

    samples = rnaseq.combine_files(samples)
    with prun.start(_wres(parallel, ["gatk", "vardict"]), samples, config,
                    dirs, "rnaseq-variation") as run_parallel:
        with profile.report("RNA-seq variant calling", dirs):
            samples = rnaseq.rnaseq_variant_calling(samples, run_parallel)

    with prun.start(_wres(parallel, ["samtools", "fastqc", "qualimap",
                                     "kraken", "gatk", "preseq"], ensure_mem={"qualimap": 4}),
                    samples, config, dirs, "qc") as run_parallel:
        with profile.report("quality control", dirs):
            samples = qcsummary.generate_parallel(samples, run_parallel)
        with profile.report("upload", dirs):
            samples = run_parallel("upload_samples", samples)
            for sample in samples:
                run_parallel("upload_samples_project", [sample])
        with profile.report("bcbioRNAseq loading", dirs):
            run_parallel("run_bcbiornaseqload", [sample])
    logger.info("Timing: finished")
    return samples

def fastrnaseqpipeline(config, run_info_yaml, parallel, dirs, samples):
    samples = rnaseq_prep_samples(config, run_info_yaml, parallel, dirs, samples)
    ww = initialize_watcher(samples)
    with prun.start(_wres(parallel, ["samtools"]), samples, config,
                    dirs, "fastrnaseq") as run_parallel:
        with profile.report("fastrnaseq", dirs):
            samples = rnaseq.fast_rnaseq(samples, run_parallel)
            ww.report("fastrnaseq", samples)
        samples = rnaseq.combine_files(samples)
        with profile.report("quality control", dirs):
            samples = qcsummary.generate_parallel(samples, run_parallel)
            ww.report("qcsummary", samples)
        with profile.report("upload", dirs):
            samples = run_parallel("upload_samples", samples)
            for samples in samples:
                run_parallel("upload_samples_project", [samples])
    logger.info("Timing: finished")
    return samples

def singlecellrnaseqpipeline(config, run_info_yaml, parallel, dirs, samples):
    samples = rnaseq_prep_samples(config, run_info_yaml, parallel, dirs, samples)
    with prun.start(_wres(parallel, ["samtools", "rapmap"]), samples, config,
                    dirs, "singlecell-rnaseq") as run_parallel:
        with profile.report("singlecell-rnaseq", dirs):
            samples = rnaseq.singlecell_rnaseq(samples, run_parallel)
        with profile.report("quality control", dirs):
            samples = qcsummary.generate_parallel(samples, run_parallel)
        with profile.report("upload", dirs):
            samples = run_parallel("upload_samples", samples)
            for samples in samples:
                run_parallel("upload_samples_project", [samples])
    logger.info("Timing: finished")
    return samples

def smallrnaseqpipeline(config, run_info_yaml, parallel, dirs, samples):
    # causes a circular import at the top level
    from bcbio.srna.group import report as srna_report

    samples = rnaseq_prep_samples(config, run_info_yaml, parallel, dirs, samples)

    with prun.start(_wres(parallel, ["aligner", "picard", "samtools"],
                          ensure_mem={"bowtie": 8, "bowtie2": 8, "star": 2}),
                    [samples[0]], config, dirs, "alignment") as run_parallel:
        with profile.report("prepare", dirs):
            samples = run_parallel("seqcluster_prepare", [samples])
        with profile.report("seqcluster alignment", dirs):
            samples = run_parallel("srna_alignment", [samples])

    with prun.start(_wres(parallel, ["aligner", "picard", "samtools"],
                            ensure_mem={"tophat": 10, "tophat2": 10, "star": 2, "hisat2": 8}),
                    samples, config, dirs, "alignment_samples",
                    multiplier=alignprep.parallel_multiplier(samples)) as run_parallel:
        with profile.report("alignment", dirs):
            samples = run_parallel("process_alignment", samples)

    with prun.start(_wres(parallel, ["picard", "miraligner"]),
                    samples, config, dirs, "annotation") as run_parallel:
        with profile.report("small RNA annotation", dirs):
            samples = run_parallel("srna_annotation", samples)

    with prun.start(_wres(parallel, ["seqcluster", "mirge"],
                          ensure_mem={"seqcluster": 8}),
                    [samples[0]], config, dirs, "cluster") as run_parallel:
        with profile.report("cluster", dirs):
            samples = run_parallel("seqcluster_cluster", [samples])

    with prun.start(_wres(parallel, ["picard", "fastqc"]),
                    samples, config, dirs, "qc") as run_parallel:
        with profile.report("quality control", dirs):
            samples = qcsummary.generate_parallel(samples, run_parallel)
        with profile.report("report", dirs):
            srna_report(samples)
        with profile.report("upload", dirs):
            samples = run_parallel("upload_samples", samples)
            for sample in samples:
                run_parallel("upload_samples_project", [sample])

    return samples

def chipseqpipeline(config, run_info_yaml, parallel, dirs, samples):
    with prun.start(_wres(parallel, ["aligner", "picard"]),
                    samples, config, dirs, "multicore",
                    multiplier=alignprep.parallel_multiplier(samples)) as run_parallel:
        with profile.report("organize samples", dirs):
            samples = run_parallel("organize_samples", [[dirs, config, run_info_yaml,
                                                            [x[0]["description"] for x in samples]]])
        with profile.report("alignment", dirs):
            samples = run_parallel("prepare_sample", samples)
            samples = run_parallel("trim_sample", samples)
            samples = run_parallel("disambiguate_split", [samples])
            samples = run_parallel("process_alignment", samples)

        with profile.report("disambiguation", dirs):
            samples = disambiguate.resolve(samples, run_parallel)
            samples = run_parallel("clean_chipseq_alignment", samples)

    with prun.start(_wres(parallel, ["peakcaller"]),
                    samples, config, dirs, "peakcalling",
                    multiplier = peaks._get_multiplier(samples)) as run_parallel:
        with profile.report("peakcalling", dirs):
            samples = peaks.peakcall_prepare(samples, run_parallel)

    with prun.start(_wres(parallel, ["picard", "fastqc"]),
                    samples, config, dirs, "qc") as run_parallel:
        with profile.report("quality control", dirs):
            samples = qcsummary.generate_parallel(samples, run_parallel)
        with profile.report("upload", dirs):
            samples = run_parallel("upload_samples", samples)
            for sample in samples:
                run_parallel("upload_samples_project", [sample])
    logger.info("Timing: finished")
    return samples


def wgbsseqpipeline(config, run_info_yaml, parallel, dirs, samples):
    with prun.start(_wres(parallel, ["fastqc", "picard"], ensure_mem={"fastqc" : 4}),
                    samples, config, dirs, "trimming") as run_parallel:
        with profile.report("organize samples", dirs):
            samples = run_parallel("organize_samples", [[dirs, config, run_info_yaml,
                                                            [x[0]["description"] for x in samples]]])
            samples = run_parallel("prepare_sample", samples)
            samples = run_parallel("trim_bs_sample", samples)

    with prun.start(_wres(parallel, ["aligner", "bismark", "picard", "samtools"]),
                    samples, config, dirs, "multicore",
                    multiplier=alignprep.parallel_multiplier(samples)) as run_parallel:
        with profile.report("alignment", dirs):
            samples = run_parallel("process_alignment", samples)

    with prun.start(_wres(parallel, ['samtools']), samples, config, dirs,
                    'deduplication') as run_parallel:
        with profile.report('deduplicate', dirs):
            samples = run_parallel('deduplicate_bismark', samples)

    with prun.start(_wres(parallel, ["caller"], ensure_mem={"caller": 5}),
                    samples, config, dirs, "multicore2",
                    multiplier=24) as run_parallel:
        with profile.report("cpg calling", dirs):
            samples = run_parallel("cpg_calling", samples)

    # with prun.start(_wres(parallel, ["picard", "fastqc", "samtools"]),
    #                 samples, config, dirs, "qc") as run_parallel:
    #     with profile.report("quality control", dirs):
    #         samples = qcsummary.generate_parallel(samples, run_parallel)
    return samples


def rnaseq_prep_samples(config, run_info_yaml, parallel, dirs, samples):
    """
    organizes RNA-seq and small-RNAseq samples, converting from BAM if
    necessary and trimming if necessary
    """
    pipeline = dd.get_in_samples(samples, dd.get_analysis)
    trim_reads_set = any([tz.get_in(["algorithm", "trim_reads"], d) for d in dd.sample_data_iterator(samples)])
    resources = ["picard"]
    needs_trimming = (_is_smallrnaseq(pipeline) or trim_reads_set)
    if needs_trimming:
        resources.append("atropos")
    with prun.start(_wres(parallel, resources),
                    samples, config, dirs, "trimming",
                    max_multicore=1 if not needs_trimming else None) as run_parallel:
        with profile.report("organize samples", dirs):
            samples = run_parallel("organize_samples", [[dirs, config, run_info_yaml,
                                                            [x[0]["description"] for x in samples]]])
            samples = run_parallel("prepare_sample", samples)
        if needs_trimming:
            with profile.report("adapter trimming", dirs):
                if _is_smallrnaseq(pipeline):
                    samples = run_parallel("trim_srna_sample", samples)
                else:
                    samples = run_parallel("trim_sample", samples)
    return samples

def _get_pipeline(item):
    from bcbio.log import logger
    analysis_type = item.get("analysis", "").lower()
    if analysis_type not in SUPPORTED_PIPELINES:
        logger.error("Cannot determine which type of analysis to run, "
                      "set in the run_info under details.")
        sys.exit(1)
    else:
        return SUPPORTED_PIPELINES[analysis_type]

def _pair_samples_with_pipelines(run_info_yaml, config):
    """Map samples defined in input file to pipelines to run.
    """
    samples = config_utils.load_config(run_info_yaml)
    if isinstance(samples, dict):
        resources = samples.pop("resources")
        samples = samples["details"]
    else:
        resources = {}
    ready_samples = []
    for sample in samples:
        if "files" in sample:
            del sample["files"]
        # add any resources to this item to recalculate global configuration
        usample = copy.deepcopy(sample)
        usample.pop("algorithm", None)
        if "resources" not in usample:
            usample["resources"] = {}
        for prog, pkvs in resources.items():
            if prog not in usample["resources"]:
                usample["resources"][prog] = {}
            if pkvs is not None:
                for key, val in pkvs.items():
                    usample["resources"][prog][key] = val
        config = config_utils.update_w_custom(config, usample)
        sample["resources"] = {}
        ready_samples.append(sample)
    paired = [(x, _get_pipeline(x)) for x in ready_samples]
    d = defaultdict(list)
    for x in paired:
        d[x[1]].append([x[0]])
    return d, config

SUPPORTED_PIPELINES = {"variant2": variant2pipeline,
                       "snp calling": variant2pipeline,
                       "variant": variant2pipeline,
                       "standard": standardpipeline,
                       "minimal": standardpipeline,
                       "rna-seq": rnaseqpipeline,
                       "smallrna-seq": smallrnaseqpipeline,
                       "chip-seq": chipseqpipeline,
                       "wgbs-seq": wgbsseqpipeline,
                       "fastrna-seq": fastrnaseqpipeline,
                       "scrna-seq": singlecellrnaseqpipeline}

def _is_smallrnaseq(pipeline):
    return pipeline.lower() == "smallrna-seq"
