"""Main entry point for distributed next-gen sequencing pipelines.

Handles running the full pipeline based on instructions
"""
import abc
from collections import defaultdict
import copy
import os
import sys
import resource
import tempfile

from bcbio import log, structural, utils, upload
from bcbio.distributed import prun
from bcbio.illumina import flowcell
from bcbio.log import logger
from bcbio.ngsalign import alignprep
from bcbio.pipeline import (archive, disambiguate, region, run_info, qcsummary,
                            rnaseq)
from bcbio.pipeline.config_utils import load_system_config
from bcbio.provenance import diagnostics, programs, profile, system, versioncheck
from bcbio.variation import coverage, ensemble, genotype, population, validate, joint

def run_main(workdir, config_file=None, fc_dir=None, run_info_yaml=None,
             parallel=None, workflow=None):
    """Run variant analysis, handling command line options.
    """
    os.chdir(workdir)
    config, config_file = load_system_config(config_file, workdir)
    if config.get("log_dir", None) is None:
        config["log_dir"] = os.path.join(workdir, "log")
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
    parallel = log.create_base_logger(config, parallel)
    log.setup_local_logging(config, parallel)
    dirs = setup_directories(work_dir, fc_dir, config, config_file)
    config_file = os.path.join(dirs["config"], os.path.basename(config_file))
    samples = run_info.organize(dirs, config, run_info_yaml)
    pipelines = _pair_samples_with_pipelines(samples)
    final = []
    with utils.curdir_tmpdir({"config": config}) as tmpdir:
        tempfile.tempdir = tmpdir
        for pipeline, pipeline_items in pipelines.items():
            pipeline_items = _add_provenance(pipeline_items, dirs, parallel, config)
            versioncheck.testall(pipeline_items)
            for xs in pipeline.run(config, config_file, parallel, dirs, pipeline_items):
                if len(xs) == 1:
                    upload.from_sample(xs[0])
                    final.append(xs[0])

def setup_directories(work_dir, fc_dir, config, config_file):
    fastq_dir, galaxy_dir, config_dir = _get_full_paths(flowcell.get_fastq_dir(fc_dir)
                                                        if fc_dir else None,
                                                        config, config_file)
    return {"fastq": fastq_dir, "galaxy": galaxy_dir,
            "work": work_dir, "flowcell": fc_dir, "config": config_dir}

def _add_provenance(items, dirs, parallel, config):
    p = programs.write_versions(dirs, config, is_wrapper=parallel.get("wrapper") is not None)
    p_db = diagnostics.initialize(dirs)
    system.write_info(dirs, parallel, config)
    out = []
    for item in items:
        entity_id = diagnostics.store_entity(item)
        item["config"]["resources"]["program_versions"] = p
        item["provenance"] = {"programs": p, "entity": entity_id,
                              "db": p_db}
        out.append([item])
    return out

# ## Utility functions


def _get_full_paths(fastq_dir, config, config_file):
    """Retrieve full paths for directories in the case of relative locations.
    """
    if fastq_dir:
        fastq_dir = utils.add_full_path(fastq_dir)
    config_dir = utils.add_full_path(os.path.dirname(config_file))
    galaxy_config_file = utils.add_full_path(config.get("galaxy_config", "universe_wsgi.ini"),
                                             config_dir)
    return fastq_dir, os.path.dirname(galaxy_config_file), config_dir

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

class AbstractPipeline:
    """
    Implement this class to participate in the Pipeline abstraction.
    name: the analysis name in the run_info.yaml file:
        design:
            - analysis: name
    run: the steps run to perform the analyses
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def name(self):
        return

    @abc.abstractmethod
    def run(self, config, config_file, parallel, dirs, samples):
        return

class Variant2Pipeline(AbstractPipeline):
    """Streamlined variant calling pipeline for large files.
    This is less generalized but faster in standard cases.
    The goal is to replace the base variant calling approach.
    """
    name = "variant2"

    @classmethod
    def run(self, config, config_file, parallel, dirs, samples):
        ## Alignment and preparation requiring the entire input file (multicore cluster)
        with prun.start(_wres(parallel, ["aligner", "samtools", "sambamba"],
                              (["reference", "fasta"], ["reference", "aligner"], ["files"])),
                        samples, config, dirs, "multicore",
                        multiplier=alignprep.parallel_multiplier(samples)) as run_parallel:
            with profile.report("alignment preparation", dirs):
                samples = run_parallel("prep_align_inputs", samples)
                samples = disambiguate.split(samples)
            with profile.report("alignment", dirs):
                samples = run_parallel("process_alignment", samples)
                samples = alignprep.merge_split_alignments(samples, run_parallel)
                samples = disambiguate.resolve(samples, run_parallel)
            with profile.report("callable regions", dirs):
                samples = run_parallel("postprocess_alignment", samples)
                samples = run_parallel("combine_sample_regions", samples)
                samples = region.clean_sample_data(samples)
            with profile.report("coverage", dirs):
                samples = coverage.summarize_samples(samples, run_parallel)

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
                                         "gemini", "samtools", "fastqc", "bamtools",
                                         "bcbio-variation-recall"]),
                        samples, config, dirs, "multicore2") as run_parallel:
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
                samples = structural.run(samples, run_parallel)
            with profile.report("population database", dirs):
                samples = population.prep_db_parallel(samples, run_parallel)
            with profile.report("quality control", dirs):
                samples = qcsummary.generate_parallel(samples, run_parallel)
            with profile.report("archive", dirs):
                samples = archive.compress(samples, run_parallel)
        logger.info("Timing: finished")
        return samples

def _debug_samples(i, samples):
    print "---", i, len(samples)
    for sample in (x[0] for x in samples):
        print "  ", sample["description"], sample.get("region"), \
            utils.get_in(sample, ("config", "algorithm", "variantcaller")), \
            utils.get_in(sample, ("config", "algorithm", "jointcaller")), \
            [x.get("variantcaller") for x in sample.get("variants", [])], \
            sample.get("work_bam")

class SNPCallingPipeline(Variant2Pipeline):
    """Back compatible: old name for variant analysis.
    """
    name = "SNP calling"

class VariantPipeline(Variant2Pipeline):
    """Back compatibility; old name
    """
    name = "variant"

class StandardPipeline(AbstractPipeline):
    """Minimal pipeline with alignment and QC.
    """
    name = "Standard"
    @classmethod
    def run(self, config, config_file, parallel, dirs, samples):
        ## Alignment and preparation requiring the entire input file (multicore cluster)
        with prun.start(_wres(parallel, ["aligner"]),
                        samples, config, dirs, "multicore") as run_parallel:
            with profile.report("alignment", dirs):
                samples = run_parallel("process_alignment", samples)
            with profile.report("callable regions", dirs):
                samples = run_parallel("postprocess_alignment", samples)
                samples = run_parallel("combine_sample_regions", samples)
                samples = region.clean_sample_data(samples)
        ## Quality control
        with prun.start(_wres(parallel, ["fastqc", "bamtools", "samtools", "qsignature"]),
                        samples, config, dirs, "multicore2") as run_parallel:
            with profile.report("quality control", dirs):
                samples = qcsummary.generate_parallel(samples, run_parallel)
        logger.info("Timing: finished")
        return samples

class MinimalPipeline(StandardPipeline):
    name = "Minimal"

class SailfishPipeline(AbstractPipeline):
    name = "sailfish"

    @classmethod
    def run(self, config, config_file, parallel, dirs, samples):
        with prun.start(_wres(parallel, ["picard", "AlienTrimmer"]),
                        samples, config, dirs, "trimming") as run_parallel:
            with profile.report("adapter trimming", dirs):
                samples = run_parallel("prepare_sample", samples)
                samples = run_parallel("trim_sample", samples)
            with prun.start(_wres(parallel, ["sailfish"]), samples, config, dirs,
                            "sailfish") as run_parallel:
                with profile.report("sailfish", dirs):
                    samples = run_parallel("run_sailfish", samples)
        return samples

class RnaseqPipeline(AbstractPipeline):
    name = "RNA-seq"

    @classmethod
    def run(self, config, config_file, parallel, dirs, samples):
        with prun.start(_wres(parallel, ["picard", "AlienTrimmer"]),
                        samples, config, dirs, "trimming") as run_parallel:
            with profile.report("adapter trimming", dirs):
                samples = run_parallel("prepare_sample", samples)
                samples = run_parallel("trim_sample", samples)
        with prun.start(_wres(parallel, ["aligner", "picard"],
                              ensure_mem={"tophat": 8, "tophat2": 8, "star": 40}),
                        samples, config, dirs, "multicore",
                        multiplier=alignprep.parallel_multiplier(samples)) as run_parallel:
            with profile.report("alignment", dirs):
                samples = disambiguate.split(samples)
                samples = run_parallel("process_alignment", samples)
        with prun.start(_wres(parallel, ["samtools", "cufflinks"]),
                        samples, config, dirs, "rnaseqcount") as run_parallel:
            with profile.report("disambiguation", dirs):
                samples = disambiguate.resolve(samples, run_parallel)
            with profile.report("transcript assembly", dirs):
                samples = rnaseq.assemble_transcripts(run_parallel, samples)
            with profile.report("estimate expression", dirs):
                samples = rnaseq.estimate_expression(samples, run_parallel)
        with prun.start(_wres(parallel, ["picard", "fastqc", "rnaseqc", "kraken"]),
                        samples, config, dirs, "persample") as run_parallel:
            with profile.report("quality control", dirs):
                samples = qcsummary.generate_parallel(samples, run_parallel)
        logger.info("Timing: finished")
        return samples

class ChipseqPipeline(AbstractPipeline):
    name = "chip-seq"

    @classmethod
    def run(self, config, config_file, parallel, dirs, samples):
        with prun.start(_wres(parallel, ["aligner", "picard"]),
                        samples, config, dirs, "multicore",
                        multiplier=alignprep.parallel_multiplier(samples)) as run_parallel:
            samples = run_parallel("prepare_sample", samples)
            samples = run_parallel("trim_sample", samples)
            samples = disambiguate.split(samples)
            samples = run_parallel("process_alignment", samples)
        with prun.start(_wres(parallel, ["picard", "fastqc"]),
                        samples, config, dirs, "persample") as run_parallel:
            samples = run_parallel("clean_chipseq_alignment", samples)
            samples = qcsummary.generate_parallel(samples, run_parallel)
        return samples

def _get_pipeline(item):
    from bcbio.log import logger
    SUPPORTED_PIPELINES = {x.name.lower(): x for x in
                           utils.itersubclasses(AbstractPipeline)}
    analysis_type = item.get("analysis", "").lower()
    if analysis_type not in SUPPORTED_PIPELINES:
        logger.error("Cannot determine which type of analysis to run, "
                      "set in the run_info under details.")
        sys.exit(1)
    else:
        return SUPPORTED_PIPELINES[analysis_type]

def _pair_samples_with_pipelines(samples):
    paired = [(x, _get_pipeline(x)) for x in samples]
    d = defaultdict(list)
    for x in paired:
        d[x[1]].append(x[0])
    return d
