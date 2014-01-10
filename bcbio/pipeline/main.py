"""Main entry point for distributed next-gen sequencing pipelines.

Handles running the full pipeline based on instructions
"""
import abc
from collections import defaultdict
import os
import sys
import argparse
import resource
import tempfile

from bcbio import install, log, structural, utils, upload
from bcbio.bam import callable
from bcbio.distributed.messaging import parallel_runner
from bcbio.distributed.ipython import global_parallel
from bcbio.log import logger
from bcbio.ngsalign import alignprep
from bcbio.pipeline import (disambiguate, lane, region, run_info, qcsummary, version,
                            rnaseq)
from bcbio.pipeline.config_utils import load_system_config
from bcbio.provenance import programs, system, versioncheck
from bcbio.server import main as server_main
from bcbio.solexa.flowcell import get_fastq_dir
from bcbio.variation.realign import parallel_realign_sample
from bcbio.variation.genotype import parallel_variantcall, combine_multiple_callers
from bcbio.variation import coverage, ensemble, population, recalibrate, validate
from bcbio.rnaseq.count import combine_count_files

def run_main(work_dir, config_file=None, fc_dir=None, run_info_yaml=None,
             numcores=None, paralleltype=None, queue=None, scheduler=None,
             upgrade=None, profile=None, workflow=None, inputs=None,
             resources="", timeout=15, retries=None):
    """Run variant analysis, handling command line options.
    """
    config, config_file = load_system_config(config_file)
    if config.get("log_dir", None) is None:
        config["log_dir"] = os.path.join(work_dir, "log")
    paralleltype, numcores = _get_cores_and_type(numcores, paralleltype, scheduler)
    parallel = {"type": paralleltype, "cores": numcores,
                "scheduler": scheduler, "queue": queue,
                "profile": profile, "module": "bcbio.distributed",
                "resources": resources, "timeout": timeout,
                "retries": retries}
    if parallel["type"] in ["local"]:
        _setup_resources()
        _run_toplevel(config, config_file, work_dir, parallel,
                      fc_dir, run_info_yaml)
    elif parallel["type"] == "ipython":
        assert parallel["queue"] is not None, "IPython parallel requires a specified queue (-q)"
        assert parallel["scheduler"] is not None, "IPython parallel requires a specified scheduler (-s)"
        _run_toplevel(config, config_file, work_dir, parallel,
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

def _get_cores_and_type(numcores, paralleltype, scheduler):
    """Return core and parallelization approach from command line providing sane defaults.
    """
    if scheduler is not None:
        paralleltype = "ipython"
    if paralleltype is None:
        paralleltype = "local"
    if numcores is None:
        numcores = 1
    return paralleltype, int(numcores)

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
    fastq_dir, galaxy_dir, config_dir = _get_full_paths(get_fastq_dir(fc_dir)
                                                        if fc_dir else None,
                                                        config, config_file)
    config_file = os.path.join(config_dir, os.path.basename(config_file))
    dirs = {"fastq": fastq_dir, "galaxy": galaxy_dir,
            "work": work_dir, "flowcell": fc_dir, "config": config_dir}
    run_items = run_info.organize(dirs, config, run_info_yaml)
    run_parallel = parallel_runner(parallel, dirs, config, config_file)

    # process each flowcell lane
    lane_items = lane.process_all_lanes(run_items, run_parallel)
    pipelines = _pair_lanes_with_pipelines(lane_items)
    final = []
    with utils.curdir_tmpdir() as tmpdir:
        tempfile.tempdir = tmpdir
        for pipeline, pipeline_items in pipelines.items():
            pipeline_items = _add_provenance(pipeline_items, dirs, run_parallel, parallel, config)
            versioncheck.testall(pipeline_items)
            for xs in pipeline.run(config, config_file, run_parallel, parallel, dirs, pipeline_items):
                if len(xs) == 1:
                    upload.from_sample(xs[0])
                    final.append(xs[0])

def _add_provenance(items, dirs, run_parallel, parallel, config):
    p = programs.write_versions(dirs, config)
    system.write_info(dirs, run_parallel, parallel, config)
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
                "server": server_main.add_subparser}
    parser = argparse.ArgumentParser(
        description="Best-practice pipelines for fully automated high throughput sequencing analysis.")
    sub_cmd = None
    if len(in_args) > 0 and in_args[0] in sub_cmds:
        subparsers = parser.add_subparsers(help="bcbio-nextgen supplemental commands")
        sub_parser = sub_cmds[in_args[0]](subparsers)
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
                                  "Supports SGE and SLURM parameters."),
                            default=[], action="append")
        parser.add_argument("--timeout", help="Number of minutes before cluster startup times out. Defaults to 15",
                            default=15, type=int)
        parser.add_argument("--retries",
                            help=("Number of retries of failed tasks during distributed processing. "
                                  "Default 0 (no retries)"),
                            default=0, type=int)
        parser.add_argument("-p", "--profile", help="Profile name to use for ipython parallel",
                            default="bcbio_nextgen")
        parser.add_argument("-w", "--workflow", help="Run a workflow with the given commandline arguments")
        parser.add_argument("-v", "--version", help="Print current version",
                            action="store_true")
    args = parser.parse_args(in_args)
    if hasattr(args, "global_config"):
        error_msg = _sanity_check_args(args)
        if error_msg:
            parser.error(error_msg)
        kwargs = {"numcores": args.numcores if args.numcores > 0 else None,
                  "paralleltype": args.paralleltype,
                  "scheduler": args.scheduler,
                  "queue": args.queue,
                  "timeout": args.timeout,
                  "retries": args.retries,
                  "resources": args.resources,
                  "profile": args.profile,
                  "workflow": args.workflow}
        kwargs = _add_inputs_to_kwargs(args, kwargs, parser)
    else:
        assert sub_cmd is not None
        if sub_cmd == "upgrade":
            # Nothing specified to do, print help
            if not args.tools and not args.tooldir and not args.install_data and args.upgrade == "skip":
                sub_parser.print_help()
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
    def run(self, config, config_file, run_parallel, parallel, dirs, lanes):
        return


class VariantPipeline(AbstractPipeline):
    name = "variant"

    @classmethod
    def run(self, config, config_file, run_parallel, parallel, dirs, lane_items):
        raise NotImplementedError("`variant` processing is deprecated: please use `variant2`"
                                  "The next version will alias variant to the new variant2 pipeline")
        lane_items = run_parallel("trim_lane", lane_items)
        align_items = run_parallel("process_alignment", lane_items)
        # process samples, potentially multiplexed across multiple lanes
        samples = organize_samples(align_items, dirs, config_file)
        samples = run_parallel("merge_sample", samples)
        samples = run_parallel("prep_recal", samples)
        samples = recalibrate.parallel_write_recal_bam(samples, run_parallel)
        samples = parallel_realign_sample(samples, run_parallel)
        samples = parallel_variantcall(samples, run_parallel)
        samples = run_parallel("postprocess_variants", samples)
        samples = combine_multiple_callers(samples)
        samples = ensemble.combine_calls_parallel(samples, run_parallel)
        samples = run_parallel("detect_sv", samples)
        samples = qcsummary.generate_parallel(samples, run_parallel)
        run_parallel("generate_bigwig", samples, {"programs": ["ucsc_bigwig"]})
        return samples

class Variant2Pipeline(AbstractPipeline):
    """Streamlined variant calling pipeline for large files.
    This is less generalized but faster in standard cases.
    The goal is to replace the base variant calling approach.
    """
    name = "variant2"

    @classmethod
    def run(self, config, config_file, run_parallel, parallel, dirs, samples):
        ## Alignment and preparation requiring the entire input file (multicore cluster)
        with global_parallel(parallel, "multicore", ["process_alignment", "postprocess_alignment"],
                             samples, dirs, config,
                             multiplier=alignprep.parallel_multiplier(samples)) as parallel:
            run_parallel = parallel_runner(parallel, dirs, config)
            logger.info("Timing: alignment")
            samples = run_parallel("prep_align_inputs", samples)
            samples = run_parallel("process_alignment", samples)
            samples = alignprep.merge_split_alignments(samples, run_parallel)
            samples = run_parallel("postprocess_alignment", samples)
            regions = callable.combine_sample_regions(samples)
            samples = region.add_region_info(samples, regions)
            samples = region.clean_sample_data(samples)
            logger.info("Timing: coverage")
            samples = coverage.summarize_samples(samples, run_parallel)

        ## Variant calling on sub-regions of the input file (full cluster)
        with global_parallel(parallel, "full", ["piped_bamprep", "variantcall_sample"],
                             samples, dirs, config,
                             multiplier=len(regions["analysis"]), max_multicore=1) as parallel:
            run_parallel = parallel_runner(parallel, dirs, config)
            logger.info("Timing: alignment post-processing")
            samples = region.parallel_prep_region(samples, regions, run_parallel)
            logger.info("Timing: variant calling")
            samples = region.parallel_variantcall_region(samples, run_parallel)

        ## Finalize variants (per-sample cluster)
        with global_parallel(parallel, "persample", ["postprocess_variants"],
                             samples, dirs, config) as parallel:
            run_parallel = parallel_runner(parallel, dirs, config)
            logger.info("Timing: variant post-processing")
            samples = run_parallel("postprocess_variants", samples)
            logger.info("Timing: validation")
            samples = run_parallel("compare_to_rm", samples)
            samples = combine_multiple_callers(samples)
            logger.info("Timing: ensemble calling")
            samples = ensemble.combine_calls_parallel(samples, run_parallel)
            samples = validate.summarize_grading(samples)
        ## Finalizing BAMs and population databases, handle multicore computation
        with global_parallel(parallel, "multicore2", ["prep_gemini_db", "delayed_bam_merge"],
                             samples, dirs, config) as parallel:
            run_parallel = parallel_runner(parallel, dirs, config)
            logger.info("Timing: prepped BAM merging")
            samples = region.delayed_bamprep_merge(samples, run_parallel)
            logger.info("Timing: structural variation")
            samples = structural.run(samples, run_parallel)
            logger.info("Timing: population database")
            samples = population.prep_db_parallel(samples, run_parallel)
            logger.info("Timing: quality control")
            samples = qcsummary.generate_parallel(samples, run_parallel)
        logger.info("Timing: finished")
        return samples

class SNPCallingPipeline(VariantPipeline):
    """Back compatible: old name for variant analysis.
    """
    name = "SNP calling"

class StandardPipeline(AbstractPipeline):
    """Minimal pipeline with alignment and QC.
    """
    name = "Standard"
    @classmethod
    def run(self, config, config_file, run_parallel, parallel, dirs, lane_items):
        ## Alignment and preparation requiring the entire input file (multicore cluster)
        with global_parallel(parallel, "multicore", ["process_alignment"],
                             lane_items, dirs, config) as parallel:
            run_parallel = parallel_runner(parallel, dirs, config)
            logger.info("Timing: alignment")
            samples = run_parallel("process_alignment", lane_items)
        ## Finalize (per-sample cluster)
        with global_parallel(parallel, "persample", ["postprocess_variants"],
                             samples, dirs, config) as parallel:
            run_parallel = parallel_runner(parallel, dirs, config)
            logger.info("Timing: quality control")
            samples = qcsummary.generate_parallel(samples, run_parallel)
        logger.info("Timing: finished")
        return samples

class MinimalPipeline(StandardPipeline):
    name = "Minimal"

class RnaseqPipeline(AbstractPipeline):
    name = "RNA-seq"

    @classmethod
    def run(self, config, config_file, run_parallel, parallel, dirs, lane_items):
        lane_items = run_parallel("trim_lane", lane_items)
        samples = disambiguate.split(lane_items)
        samples = run_parallel("process_alignment", samples)
        samples = disambiguate.resolve(samples, run_parallel)
        samples = rnaseq.estimate_expression(samples, run_parallel)
        combined = combine_count_files([x[0].get("count_file") for x in samples])
        for x in samples:
            x[0]["combined_counts"] = combined

        samples = qcsummary.generate_parallel(samples, run_parallel)
        #run_parallel("generate_bigwig", samples, {"programs": ["ucsc_bigwig"]})
        return samples

class ChipseqPipeline(AbstractPipeline):
    name = "chip-seq"

    @classmethod
    def run(self, config, config_file, run_parallel, parallel, dirs, lane_items):
        lane_items = run_parallel("trim_lane", lane_items)
        samples = disambiguate.split(lane_items)
        samples = run_parallel("process_alignment", samples)
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
