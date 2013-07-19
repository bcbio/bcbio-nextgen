"""Main entry point for distributed next-gen sequencing pipelines.

Handles running the full pipeline based on instructions
"""
import abc
import os
import sys
import argparse
from collections import defaultdict


from bcbio import install, log, utils, upload
from bcbio.bam import callable
from bcbio import rnaseq
from bcbio.distributed.messaging import parallel_runner
from bcbio.distributed.ipython import global_parallel
from bcbio.log import logger
from bcbio.pipeline.run_info import get_run_info
from bcbio.pipeline.demultiplex import add_multiplex_across_lanes
from bcbio.pipeline.merge import organize_samples
from bcbio.pipeline import lane, region, qcsummary
from bcbio.provenance import programs
from bcbio.solexa.flowcell import get_fastq_dir
from bcbio.variation.realign import parallel_realign_sample
from bcbio.variation.genotype import parallel_variantcall, combine_multiple_callers
from bcbio.variation import ensemble, population, recalibrate, validate

def run_main(config, config_file, work_dir, parallel,
         fc_dir=None, run_info_yaml=None):
    """
    Run toplevel analysis, processing a set of input files.
    config_file -- Main YAML configuration file with system parameters
    fc_dir -- Directory of fastq files to process
    run_info_yaml -- YAML configuration file specifying inputs to process
    """
    parallel = log.create_base_logger(config, parallel)
    log.setup_local_logging(config, parallel)
    fc_name, fc_date, run_info = get_run_info(fc_dir, config, run_info_yaml)
    fastq_dir, galaxy_dir, config_dir = _get_full_paths(get_fastq_dir(fc_dir)
                                                        if fc_dir else None,
                                                        config, config_file)
    config_file = os.path.join(config_dir, os.path.basename(config_file))
    dirs = {"fastq": fastq_dir, "galaxy": galaxy_dir,
            "work": work_dir, "flowcell": fc_dir, "config": config_dir}
    run_parallel = parallel_runner(parallel, dirs, config, config_file)

    # process each flowcell lane
    run_items = add_multiplex_across_lanes(run_info["details"],
                                           dirs["fastq"], fc_name)
    lanes = ((info, fc_name, fc_date, dirs, config) for info in run_items)
    lane_items = lane.process_all_lanes(lanes, run_parallel)
    pipelines = _pair_lanes_with_pipelines(lane_items)
    for pipeline, pipeline_items in pipelines.items():
        pipeline_items = _add_provenance(pipeline_items, dirs, config)
        for xs in pipeline.run(config, config_file, run_parallel, parallel, dirs, pipeline_items):
            if len(xs) == 1:
                upload.from_sample(xs[0])
    qcsummary.write_metrics(run_info, fc_name, fc_date, dirs)

def _add_provenance(items, dirs, config):
    p = programs.write_versions(dirs, config)
    out = []
    for item_info in items:
        item = item_info[2]
        if item.get("upload"):
            entity_id = "%s.%s.%s" % (item["upload"]["fc_date"],
                                      item["upload"]["fc_name"],
                                      item["description"])
        else:
            entity_id = item["description"]
        item["provenance"] = {"programs": p, "entity": entity_id}
        cur = list(item_info)
        cur[2] = item
        out.append(cur)
    return out

# ## Utility functions

def parse_cl_args(in_args):
    """Parse input commandline arguments, handling multiple cases.

    Returns the main config file and set of kwargs.
    """
    parser = argparse.ArgumentParser(
        description= "Best-practice pipelines for fully automated high throughput sequencing analysis.")
    if in_args[0] in ["upgrade"]:
        subparsers = parser.add_subparsers(help="bcbio-nextgen supplemental commands")
        install.add_subparser(subparsers)
    else:
        parser.add_argument("global_config", help="Global YAML configuration file specifying details about the system",
                            nargs="?")
        parser.add_argument("fc_dir", help="A directory of Illumina output or fastq files to process (optional)",
                            nargs="?")
        parser.add_argument("run_config", help="YAML file with details about samples to process (optional)",
                            nargs="*")
        parser.add_argument("-n", "--numcores", type=int, default=0)
        parser.add_argument("-t", "--paralleltype", help="Approach to parallelization",
                            choices=["local", "ipython", "messaging"], default="local")
        parser.add_argument("-s", "--scheduler", help="Schedulerto use for ipython parallel",
                            choices=["lsf", "sge", "torque"])
        parser.add_argument("-q", "--queue", help="Scheduler queue to run jobs on, for ipython parallel")
        parser.add_argument("-r", "--resources",
                            help=("Cluster specific resources specifications. Provide multiple specs separated by ';'\n"
                                  "Supports SGE: translated to -l parameters"),
                            default="")
        parser.add_argument("--timeout", help="Number of minutes before cluster startup times out. Defaults to 15",
                            default=15, type=int)
        parser.add_argument("--retries",
                            help=("Number of retries of failed tasks during distributed processing. "
                                  "Default 0 (no retries)"),
                            default=0, type=int)
        parser.add_argument("-p", "--profile", help="Profile name to use for ipython parallel",
                            default="bcbio_nextgen")
        parser.add_argument("-u", "--upgrade", help="Perform an upgrade of bcbio_nextgen in place.",
                            choices = ["stable", "development", "system"])
        parser.add_argument("-w", "--workflow", help="Run a workflow with the given commandline arguments")
    args = parser.parse_args(in_args)
    if hasattr(args, "fc_dir"):
        inputs = [x for x in [args.global_config, args.fc_dir] + args.run_config
                  if x is not None]
        kwargs = {"numcores": args.numcores if args.numcores > 0 else None,
                  "paralleltype": args.paralleltype,
                  "scheduler": args.scheduler,
                  "queue": args.queue,
                  "timeout": args.timeout,
                  "retries": args.retries,
                  "resources": args.resources,
                  "profile": args.profile,
                  "upgrade": args.upgrade,
                  "workflow": args.workflow,
                  "inputs": inputs}
    else:
        args.global_config = None
        args.fc_dir = None
        kwargs = {"args": args,
                  "upgrade": args.upgrade}
    if args.fc_dir is not None and len(args.run_config) == 1:
        kwargs["fc_dir"] = args.fc_dir
        kwargs["run_info_yaml"] = args.run_config[0]
    elif args.fc_dir is not None:
        if os.path.isfile(args.fc_dir):
            kwargs["run_info_yaml"] = args.fc_dir
        else:
            kwargs["fc_dir"] = args.fc_dir
    elif args.upgrade is None and args.workflow is None:
        parser.print_help()
        sys.exit()
    return args.global_config, kwargs

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
    def run(self, config, config_file, run_parallel, parallel, dirs, lane_items):
        ## Alignment and preparation requiring the entire input file (multicore cluster)
        with global_parallel(parallel, "multicore", ["align_prep_full"],
                             lane_items, dirs["work"], config) as parallel:
            run_parallel = parallel_runner(parallel, dirs, config)
            logger.info("Timing: alignment")
            samples = run_parallel("align_prep_full", [list(x) + [config_file] for x in lane_items])
            regions = callable.combine_sample_regions(samples)
            samples = region.add_region_info(samples, regions)
            samples = region.clean_sample_data(samples)

        ## Variant calling on sub-regions of the input file (full cluster)
        with global_parallel(parallel, "full", ["piped_bamprep", "variantcall_sample"],
                             samples, dirs["work"], config,
                             multiplier=len(regions["analysis"])) as parallel:
            run_parallel = parallel_runner(parallel, dirs, config)
            logger.info("Timing: alignment post-processing")
            samples = region.parallel_prep_region(samples, regions, run_parallel)
            logger.info("Timing: variant calling")
            samples = region.parallel_variantcall_region(samples, run_parallel)

        ## Finalize variants (per-sample cluster)
        with global_parallel(parallel, "persample", ["postprocess_variants"],
                             samples, dirs["work"], config) as parallel:
            run_parallel = parallel_runner(parallel, dirs, config)
            logger.info("Timing: variant post-processing")
            samples = run_parallel("postprocess_variants", samples)
            samples = combine_multiple_callers(samples)
            logger.info("Timing: ensemble calling")
            samples = ensemble.combine_calls_parallel(samples, run_parallel)
            logger.info("Timing: prepped BAM merging")
            samples = region.delayed_bamprep_merge(samples, run_parallel)
            logger.info("Timing: validation")
            samples = run_parallel("compare_to_rm", samples)
            samples = validate.summarize_grading(samples)
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

class MinimalPipeline(VariantPipeline):
    name = "Minimal"

class StandardPipeline(VariantPipeline):
    name = "Standard"

class RnaseqPipeline(AbstractPipeline):
    name = "RNA-seq"

    @classmethod
    def run(self, config, config_file, run_parallel, parallel, dirs, lane_items):
        lane_items = run_parallel("trim_lane", lane_items)
        align_items = run_parallel("process_alignment", lane_items)
        # process samples, potentially multiplexed across multiple lanes
        samples = organize_samples(align_items, dirs, config_file)
        samples = run_parallel("merge_sample", samples)
        samples = run_parallel("generate_transcript_counts", samples)
        samples = qcsummary.generate_parallel(samples, run_parallel)
        samples = rnaseq.qc.sample_summary(samples)
        #run_parallel("generate_bigwig", samples, {"programs": ["ucsc_bigwig"]})
        return samples

def _get_pipeline(lane_item):
    from bcbio.log import logger
    SUPPORTED_PIPELINES = {x.name: x for x in
                           utils.itersubclasses(AbstractPipeline)}
    analysis_type = lane_item[2].get("analysis")
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
