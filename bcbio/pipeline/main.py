"""Main entry point for distributed next-gen sequencing pipelines.

Handles running the full pipeline based on instructions
"""
import abc
from collections import defaultdict
import copy
import os
import sys
import argparse
import resource
import tempfile

from bcbio import install, log, structural, utils, upload
from bcbio.distributed import clargs, prun, runfn
from bcbio.illumina import flowcell, machine
from bcbio.log import logger
from bcbio.ngsalign import alignprep
from bcbio.pipeline import (disambiguate, region, run_info, qcsummary,
                            version, rnaseq)
from bcbio.pipeline.config_utils import load_system_config
from bcbio.provenance import programs, profile, system, versioncheck
from bcbio.server import main as server_main
from bcbio.variation.genotype import combine_multiple_callers
from bcbio.variation import coverage, ensemble, population, validate
from bcbio.rnaseq.count import (combine_count_files,
                                annotate_combined_count_file)

def run_main(workdir, config_file=None, fc_dir=None, run_info_yaml=None,
             parallel=None, workflow=None):
    """Run variant analysis, handling command line options.
    """
    os.chdir(workdir)
    config, config_file = load_system_config(config_file, workdir)
    if config.get("log_dir", None) is None:
        config["log_dir"] = os.path.join(workdir, "log")
    if parallel["type"] in ["local"]:
        _setup_resources()
        _run_toplevel(config, config_file, workdir, parallel,
                      fc_dir, run_info_yaml)
    elif parallel["type"] == "ipython":
        assert parallel["queue"] is not None, "IPython parallel requires a specified queue (-q)"
        assert parallel["scheduler"] is not None, "IPython parallel requires a specified scheduler (-s)"
        _run_toplevel(config, config_file, workdir, parallel,
                      fc_dir, run_info_yaml)
    else:
        raise ValueError("Unexpected type of parallel run: %s" % parallel["type"])

def _setup_resources():
    """Attempt to increase resource limits up to hard limits.

    This allows us to avoid out of file handle limits where we can
    move beyond the soft limit up to the hard limit.
    """
    target_procs = 50000
    cur_proc, max_proc = resource.getrlimit(resource.RLIMIT_NPROC)
    target_proc = min(max_proc, target_procs)
    resource.setrlimit(resource.RLIMIT_NPROC, (max(cur_proc, target_proc), max_proc))
    cur_hdls, max_hdls = resource.getrlimit(resource.RLIMIT_NOFILE)
    target_hdls = min(max_hdls, target_procs)
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
    fastq_dir, galaxy_dir, config_dir = _get_full_paths(flowcell.get_fastq_dir(fc_dir)
                                                        if fc_dir else None,
                                                        config, config_file)
    config_file = os.path.join(config_dir, os.path.basename(config_file))
    dirs = {"fastq": fastq_dir, "galaxy": galaxy_dir,
            "work": work_dir, "flowcell": fc_dir, "config": config_dir}
    samples = run_info.organize(dirs, config, run_info_yaml)
    pipelines = _pair_lanes_with_pipelines(samples)
    final = []
    with utils.curdir_tmpdir() as tmpdir:
        tempfile.tempdir = tmpdir
        for pipeline, pipeline_items in pipelines.items():
            pipeline_items = _add_provenance(pipeline_items, dirs, parallel, config)
            versioncheck.testall(pipeline_items)
            for xs in pipeline.run(config, config_file, parallel, dirs, pipeline_items):
                if len(xs) == 1:
                    upload.from_sample(xs[0])
                    final.append(xs[0])

def _add_provenance(items, dirs, parallel, config):
    p = programs.write_versions(dirs, config, is_wrapper=parallel.get("wrapper") is not None)
    system.write_info(dirs, parallel, config)
    out = []
    for item in items:
        if item.get("upload") and item["upload"].get("fc_name"):
            entity_id = "%s.%s.%s" % (item["upload"]["fc_date"],
                                      item["upload"]["fc_name"],
                                      item["description"])
        else:
            entity_id = item["description"]
        item["config"]["resources"]["program_versions"] = p
        item["provenance"] = {"programs": p, "entity": entity_id}
        out.append([item])
    return out

# ## Utility functions

def _sanity_check_args(args):
    """Ensure dependent arguments are correctly specified
    """
    if "scheduler" in args and "queue" in args:
        if args.scheduler and not args.queue:
            return "IPython parallel scheduler (-s) specified. This also requires a queue (-q)."
        elif args.queue and not args.scheduler:
            return "IPython parallel queue (-q) supplied. This also requires a scheduler (-s)."
        elif args.paralleltype == "ipython" and (not args.queue or not args.scheduler):
            return "IPython parallel requires queue (-q) and scheduler (-s) arguments."

def parse_cl_args(in_args):
    """Parse input commandline arguments, handling multiple cases.

    Returns the main config file and set of kwargs.
    """
    sub_cmds = {"upgrade": install.add_subparser,
                "server": server_main.add_subparser,
                "runfn": runfn.add_subparser,
                "version": programs.add_subparser,
                "sequencer": machine.add_subparser}
    parser = argparse.ArgumentParser(
        description="Best-practice pipelines for fully automated high throughput sequencing analysis.")
    sub_cmd = None
    if len(in_args) > 0 and in_args[0] in sub_cmds:
        subparsers = parser.add_subparsers(help="bcbio-nextgen supplemental commands")
        sub_cmds[in_args[0]](subparsers)
        sub_cmd = in_args[0]
    else:
        parser.add_argument("global_config", help="Global YAML configuration file specifying details "
                            "about the system (optional, defaults to installed bcbio_system.yaml)",
                            nargs="?")
        parser.add_argument("fc_dir", help="A directory of Illumina output or fastq files to process (optional)",
                            nargs="?")
        parser.add_argument("run_config", help="YAML file with details about samples to process "
                            "(required, unless using Galaxy LIMS as input)",
                            nargs="*")
        parser.add_argument("-n", "--numcores", help="Total cores to use for processing",
                            type=int, default=1)
        parser.add_argument("-t", "--paralleltype", help="Approach to parallelization",
                            choices=["local", "ipython"], default="local")
        parser.add_argument("-s", "--scheduler", help="Scheduler to use for ipython parallel",
                            choices=["lsf", "sge", "torque", "slurm"])
        parser.add_argument("-q", "--queue", help="Scheduler queue to run jobs on, for ipython parallel")
        parser.add_argument("-r", "--resources",
                            help=("Cluster specific resources specifications. Can be specified multiple times.\n"
                                  "Supports SGE, Torque, LSF and SLURM parameters."),
                            default=[], action="append")
        parser.add_argument("--timeout", help="Number of minutes before cluster startup times out. Defaults to 15",
                            default=15, type=int)
        parser.add_argument("--retries",
                            help=("Number of retries of failed tasks during distributed processing. "
                                  "Default 0 (no retries)"),
                            default=0, type=int)
        parser.add_argument("-p", "--tag", help="Tag name to label jobs on the cluster",
                            default="")
        parser.add_argument("-w", "--workflow", help="Run a workflow with the given commandline arguments")
        parser.add_argument("--workdir", help="Directory to process in. Defaults to current working directory",
                            default=os.getcwd())
        parser.add_argument("-v", "--version", help="Print current version",
                            action="store_true")
    args = parser.parse_args(in_args)
    if hasattr(args, "global_config"):
        error_msg = _sanity_check_args(args)
        if error_msg:
            parser.error(error_msg)
        kwargs = {"parallel": clargs.to_parallel(args),
                  "workflow": args.workflow,
                  "workdir": args.workdir}
        kwargs = _add_inputs_to_kwargs(args, kwargs, parser)
    else:
        assert sub_cmd is not None
        kwargs = {"args": args,
                  "config_file": None,
                  sub_cmd: True}
    return kwargs

def _add_inputs_to_kwargs(args, kwargs, parser):
    """Convert input system config, flow cell directory and sample yaml to kwargs.

    Handles back compatibility with previous commandlines while allowing flexible
    specification of input parameters.
    """
    inputs = [x for x in [args.global_config, args.fc_dir] + args.run_config
              if x is not None]
    global_config = "bcbio_system.yaml"  # default configuration if not specified
    if len(inputs) == 1:
        if os.path.isfile(inputs[0]):
            fc_dir = None
            run_info_yaml = inputs[0]
        else:
            fc_dir = inputs[0]
            run_info_yaml = None
    elif len(inputs) == 2:
        if os.path.isfile(inputs[0]):
            global_config = inputs[0]
            if os.path.isfile(inputs[1]):
                fc_dir = None
                run_info_yaml = inputs[1]
            else:
                fc_dir = inputs[1]
                run_info_yaml = None
        else:
            fc_dir, run_info_yaml = inputs
    elif len(inputs) == 3:
        global_config, fc_dir, run_info_yaml = inputs
    elif kwargs.get("workflow", "") == "template":
        kwargs["inputs"] = inputs
        return kwargs
    elif args.version:
        print version.__version__
        sys.exit()
    else:
        print "Incorrect input arguments", inputs
        parser.print_help()
        sys.exit()
    if fc_dir:
        fc_dir = os.path.abspath(fc_dir)
    if run_info_yaml:
        run_info_yaml = os.path.abspath(run_info_yaml)
    if kwargs.get("workflow"):
        kwargs["inputs"] = inputs
    kwargs["config_file"] = global_config
    kwargs["fc_dir"] = fc_dir
    kwargs["run_info_yaml"] = run_info_yaml
    return kwargs

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
    def run(self, config, config_file, parallel, dirs, lanes):
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
        with prun.start(_wres(parallel, ["aligner", "gatk"],
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
                regions = run_parallel("combine_sample_regions", [samples])[0]
                samples = region.add_region_info(samples, regions)
                samples = region.clean_sample_data(samples)
            with profile.report("coverage", dirs):
                samples = coverage.summarize_samples(samples, run_parallel)

        ## Variant calling on sub-regions of the input file (full cluster)
        with prun.start(_wres(parallel, ["gatk", "picard", "variantcaller"]),
                        samples, config, dirs, "full",
                        multiplier=len(regions["analysis"]), max_multicore=1) as run_parallel:
            with profile.report("alignment post-processing", dirs):
                samples = region.parallel_prep_region(samples, regions, run_parallel)
            with profile.report("variant calling", dirs):
                samples = region.parallel_variantcall_region(samples, run_parallel)

        ## Finalize variants (per-sample cluster)
        with prun.start(_wres(parallel, ["gatk", "gatk-vqsr", "snpeff", "bcbio_variation"]),
                        samples, config, dirs, "persample") as run_parallel:
            with profile.report("variant post-processing", dirs):
                samples = run_parallel("postprocess_variants", samples)
            with profile.report("validation", dirs):
                samples = run_parallel("compare_to_rm", samples)
                samples = combine_multiple_callers(samples)
        ## Finalizing BAMs and population databases, handle multicore computation
        with prun.start(_wres(parallel, ["gemini", "samtools", "fastqc", "bamtools", "bcbio_variation",
                                         "bcbio-variation-recall"]),
                        samples, config, dirs, "multicore2") as run_parallel:
            with profile.report("prepped BAM merging", dirs):
                samples = region.delayed_bamprep_merge(samples, run_parallel)
            with profile.report("ensemble calling", dirs):
                samples = ensemble.combine_calls_parallel(samples, run_parallel)
                samples = validate.summarize_grading(samples)
            with profile.report("structural variation", dirs):
                samples = structural.run(samples, run_parallel)
            with profile.report("population database", dirs):
                samples = population.prep_db_parallel(samples, run_parallel)
            with profile.report("quality control", dirs):
                samples = qcsummary.generate_parallel(samples, run_parallel)
        logger.info("Timing: finished")
        return samples

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
    def run(self, config, config_file, parallel, dirs, lane_items):
        ## Alignment and preparation requiring the entire input file (multicore cluster)
        with prun.start(_wres(parallel, ["aligner"]),
                        lane_items, config, dirs, "multicore") as run_parallel:
            logger.info("Timing: alignment")
            samples = run_parallel("process_alignment", lane_items)
        ## Finalize (per-sample cluster)
        with prun.start(_wres(parallel, ["fastqc", "bamtools"]),
                        samples, config, dirs, "persample") as run_parallel:
            logger.info("Timing: quality control")
            samples = qcsummary.generate_parallel(samples, run_parallel)
        logger.info("Timing: finished")
        return samples

class MinimalPipeline(StandardPipeline):
    name = "Minimal"

class RnaseqPipeline(AbstractPipeline):
    name = "RNA-seq"

    @classmethod
    def run(self, config, config_file, parallel, dirs, samples):
        with prun.start(_wres(parallel, ["picard", "AlienTrimmer"]),
                        samples, config, dirs, "trimming") as run_parallel:
            with profile.report("adapter trimming", dirs):
                samples = run_parallel("process_lane", samples)
                samples = run_parallel("trim_lane", samples)
        with prun.start(_wres(parallel, ["aligner"],
                              ensure_mem={"tophat": 8, "tophat2": 8, "star": 30}),
                        samples, config, dirs, "multicore",
                        multiplier=alignprep.parallel_multiplier(samples)) as run_parallel:
            with profile.report("alignment", dirs):
                samples = disambiguate.split(samples)
                samples = run_parallel("process_alignment", samples)

        with prun.start(_wres(parallel, ["samtools", "cufflinks"]),
                        samples, config, dirs, "rnaseqcount") as run_parallel:
            with profile.report("disambiguation", dirs):
                samples = disambiguate.resolve(samples, run_parallel)
            with profile.report("estimate expression", dirs):
                samples = rnaseq.estimate_expression(samples, run_parallel)

        combined = combine_count_files([x[0].get("count_file") for x in samples])
        gtf_file = utils.get_in(samples[0][0], ('genome_resources', 'rnaseq',
                                                'transcripts'), None)
        annotated = annotate_combined_count_file(combined, gtf_file)
        for x in samples:
            x[0]["combined_counts"] = combined
            if annotated:
                x[0]["annotated_combined_counts"] = annotated

        with prun.start(_wres(parallel, ["picard", "fastqc", "rnaseqc"]),
                        samples, config, dirs, "persample") as run_parallel:
            with profile.report("quality control", dirs):
                samples = qcsummary.generate_parallel(samples, run_parallel)
        return samples

class ChipseqPipeline(AbstractPipeline):
    name = "chip-seq"

    @classmethod
    def run(self, config, config_file, parallel, dirs, samples):
        with prun.start(_wres(parallel, ["aligner", "picard"]),
                        samples, config, dirs, "multicore",
                        multiplier=alignprep.parallel_multiplier(samples)) as run_parallel:
            samples = run_parallel("process_lane", samples)
            samples = run_parallel("trim_lane", samples)
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

def _pair_lanes_with_pipelines(lane_items):
    paired = [(x, _get_pipeline(x)) for x in lane_items]
    d = defaultdict(list)
    for x in paired:
        d[x[1]].append(x[0])
    return d
