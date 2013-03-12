"""Main entry point for distributed next-gen sequencing pipelines.

Handles running the full pipeline based on instructions
"""
import abc
import os
import sys
import argparse
from collections import defaultdict

from bcbio import utils, upload
from bcbio.bam import callable
from bcbio.log import setup_logging
from bcbio.distributed.messaging import parallel_runner
from bcbio.pipeline.run_info import get_run_info
from bcbio.pipeline.demultiplex import add_multiplex_across_lanes
from bcbio.pipeline.merge import organize_samples
from bcbio.pipeline.qcsummary import write_metrics, write_project_summary
from bcbio.pipeline import region
from bcbio.solexa.flowcell import get_fastq_dir
from bcbio.variation.realign import parallel_realign_sample
from bcbio.variation.genotype import parallel_variantcall, combine_multiple_callers
from bcbio.variation import recalibrate

def run_main(config, config_file, work_dir, parallel,
         fc_dir=None, run_info_yaml=None):
    """
    Run toplevel analysis, processing a set of input files.
    config_file -- Main YAML configuration file with system parameters
    fc_dir -- Directory of fastq files to process
    run_info_yaml -- YAML configuration file specifying inputs to process
    """

    setup_logging(config)
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
    lane_items = run_parallel("process_lane", lanes)
    pipelines = _pair_lanes_with_pipelines(lane_items)
    for pipeline, pipeline_items in pipelines.items():
        for xs in pipeline.run(config, config_file, run_parallel, dirs, pipeline_items):
            assert len(xs) == 1
            upload.from_sample(xs[0])
    write_metrics(run_info, fc_name, fc_date, dirs)

# ## Utility functions

def parse_cl_args(in_args):
    """Parse input commandline arguments, handling multiple cases.

    Returns the main config file and set of kwargs.
    """
    parser = argparse.ArgumentParser(
        description= "Best-practice pipelines for fully automated high throughput sequencing analysis.")
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
                        choices=["lsf", "sge"])
    parser.add_argument("-q", "--queue", help="Scheduler queue to run jobs on, for ipython parallel")
    parser.add_argument("-p", "--profile", help="Profile name to use for ipython parallel",
                        default="bcbio_nextgen")
    parser.add_argument("-u", "--upgrade", help="Perform an upgrade of bcbio_nextgen in place.",
                        choices = ["stable", "development", "system"])
    parser.add_argument("-w", "--workflow", help="Run a workflow with the given commandline arguments")

    args = parser.parse_args(in_args)
    inputs = [x for x in [args.global_config, args.fc_dir] + args.run_config
              if x is not None]
    kwargs = {"numcores": args.numcores if args.numcores > 0 else None,
              "paralleltype": args.paralleltype,
              "scheduler": args.scheduler,
              "queue": args.queue,
              "profile": args.profile,
              "upgrade": args.upgrade,
              "workflow": args.workflow,
              "inputs": inputs}
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
    def run(self, config, config_file, run_parallel, dirs, lanes):
        return


class VariantPipeline(AbstractPipeline):
    name = "variant"

    @classmethod
    def run(self, config, config_file, run_parallel, dirs, lane_items):
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
        samples = run_parallel("detect_sv", samples)
        samples = run_parallel("combine_calls", samples)
        samples = run_parallel("process_sample", samples)
        run_parallel("generate_bigwig", samples, {"programs": ["ucsc_bigwig"]})
        write_project_summary(samples)
        return samples

class Variant2Pipeline(AbstractPipeline):
    """Streamlined variant calling pipeline for large files.
    This is less generalized but faster in standard cases.
    The goal is to replace the base variant calling approach.
    """
    name = "variant2"

    @classmethod
    def run(self, config, config_file, run_parallel, dirs, lane_items):
        # Handle alignment and preparation requiring the entire input file
        samples = run_parallel("align_prep_full", (list(x) + [config_file] for x in lane_items))
        samples = callable.combine_sample_regions(samples)
        # Handle all variant calling on sub-regions of the input file
        samples = region.parallel_prep_region(samples, run_parallel)
        samples = region.parallel_variantcall_region(samples, run_parallel)
        samples = combine_multiple_callers(samples)
        samples = run_parallel("combine_calls", samples)
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
    def run(self, config, config_file, run_parallel, dirs, lane_items):
        lane_items = run_parallel("trim_lane", lane_items)
        align_items = run_parallel("process_alignment", lane_items)
        # process samples, potentially multiplexed across multiple lanes
        samples = organize_samples(align_items, dirs, config_file)
        samples = run_parallel("merge_sample", samples)
        samples = run_parallel("generate_transcript_counts", samples)
        run_parallel("generate_bigwig", samples, {"programs": ["ucsc_bigwig"]})
        write_project_summary(samples)
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
