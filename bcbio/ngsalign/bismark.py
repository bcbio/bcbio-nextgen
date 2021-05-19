import os
import glob
import sys
import shutil
import pysam

from bcbio.pipeline import config_utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.utils import (safe_makedir, file_exists)
from bcbio.provenance import do
from bcbio import utils
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio import bam
from bcbio import broad
from bcbio.wgbsseq import kits


def align(fastq_file, pair_file, ref_file, names, align_dir, data):
    assert data["analysis"].lower().startswith("wgbs-seq"), "No comparible alignment."
    config = data["config"]
    sample = dd.get_sample_name(data)
    out_prefix = os.path.join(align_dir, dd.get_lane(data))
    out_dir = os.path.join(align_dir, "%s_bismark" % dd.get_lane(data))

    if not ref_file:
        logger.error("bismark index not found. You can install "
                     "the index for your genome with: bcbio_nextgen.py upgrade "
                     "--aligners bismark --genomes genome-build-name --data")
        sys.exit(1)

    final_out = os.path.join(align_dir, "{0}.bam".format(sample))
    if file_exists(final_out):
        data = dd.set_work_bam(data, final_out)
        data["bam_report"] = glob.glob(os.path.join(out_dir, "*report.txt"))[0]
        data = dd.update_summary_qc(data, "bismark", base=data["bam_report"])
        return data

    bismark = config_utils.get_program("bismark", config)
    # bismark uses 5 threads/sample and ~12GB RAM/sample (hg38)
    resources = config_utils.get_resources("bismark", data["config"])
    max_cores = dd.get_num_cores(data)
    max_mem = config_utils.convert_to_bytes(resources.get("memory", "1G")) / (1024.0 * 1024.0)
    instances = calculate_bismark_instances(max_cores, max_mem * max_cores)
    # override instances if specified in the config
    if resources and resources.get("bismark_threads"):
        instances = resources.get("bismark_threads")
        logger.info(f"Using {instances} bismark instances - overriden by resources")
    bowtie_threads = 2
    if resources and resources.get("bowtie_threads"):
        bowtie_threads = resources.get("bowtie_threads")
    logger.info(f"Using {bowtie_threads} bowtie threads per bismark instance")
    kit = kits.KITS.get(dd.get_kit(data), None)
    directional = "--non_directional" if kit and not kit.is_directional else ""

    other_opts = ""
    if resources and resources.get("options", []):
        other_opts = resources.get("options", [])
        if "--non_directional" in other_opts:
            if directional != "--non_directional":
                directional = "--non_directional"
                logger.info(f"Directional setting in the kit != setting in the yaml/resources, using {directional}")
            other_opts.remove("--non_directional")
        other_opts = " ".join([str(x) for x in other_opts]).strip()

    fastq_files = " ".join([fastq_file, pair_file]) if pair_file else fastq_file
    safe_makedir(align_dir)
    cmd = "{bismark} {other_opts} {directional} --bowtie2 --temp_dir {tx_out_dir} --gzip --parallel {instances} -p {bowtie_threads} -o {tx_out_dir} --unmapped {ref_file} {fastq_file} "
    if pair_file:
        fastq_file = "-1 %s -2 %s" % (fastq_file, pair_file)
    raw_bam = glob.glob(out_dir + "/*bismark*bt2*bam")
    if not raw_bam:
        with tx_tmpdir() as tx_out_dir:
            run_message = "Running Bismark aligner on %s and %s" % (fastq_file, ref_file)
            do.run(cmd.format(**locals()), run_message, None)
            shutil.move(tx_out_dir, out_dir)
        raw_bam = glob.glob(out_dir + "/*bismark*bt2*bam")
    # don't process bam in the bismark pipeline!
    utils.symlink_plus(raw_bam[0], final_out)
    data = dd.set_work_bam(data, final_out)
    data["bam_report"] = glob.glob(os.path.join(out_dir, "*report.txt"))[0]
    data = dd.update_summary_qc(data, "bismark", base=data["bam_report"])
    return data


def _process_bam(bam_file, in_fastq, sample, reference, config):
        broad_runner = broad.runner_from_config(config)
        names = {'rg': in_fastq, 'library': 'WGBS_LIB', 'pl': 'Illumina', 'pu': 'R1', 'sm': in_fastq, 'sample': sample}
        out_fix_bam = broad_runner.run_fn("picard_fix_rgs", bam_file, names)
        order_bam = utils.append_stem(out_fix_bam, "_order")
        broad_runner.run_fn("picard_reorder", out_fix_bam, reference, order_bam)
        bam.index(order_bam, config)
        # order_bam = _set_quality(order_bam)
        # bam.index(order_bam, config)
        return order_bam

def remap_index_fn(ref_file):
    """Map sequence references to equivalent bismark indexes
    """
    return os.path.join(os.path.dirname(os.path.dirname(ref_file)), "bismark")

def _set_quality(in_bam):
    """
    change all quality to 255
    """
    bam = pysam.AlignmentFile(in_bam, "rb")
    out_file = utils.append_stem(in_bam, "_normqual")
    if file_exists(out_file):
        return out_file
    with file_transaction(out_file) as tx_out:
        with pysam.AlignmentFile(tx_out, "wb", template=bam) as out_handle:
            for read in bam.fetch():
                read.mapping_quality = 255
                out_handle.write(read)
    return out_file

def index(ref_file, out_dir, data):
    """Create a bismark index in the defined reference directory.
    """
    (ref_dir, local_file) = os.path.split(ref_file)
    gtf_file = dd.get_transcriptome_gtf(data, default=dd.get_gtf_file(data))
    bismark = config_utils.find_program("bismark", data["config"])
    if not utils.file_exists(gtf_file):
        raise ValueError("%s not found, could not create a bismark index." % (gtf_file))
    if not utils.file_exists(out_dir):
        with tx_tmpdir(data, os.path.dirname(out_dir)) as tx_out_dir:
            num_cores = dd.get_cores(data)
            other_opts = config_utils.get_resources("bismark", data["config"]).get("options", [])
            other_opts = " ".join([str(x) for x in other_opts]).strip()
            cmd = "{bismark} {other_opts} --bowtie2 -p {num_cores} -n 1 -o {tx_out_dir} --basename {sample} --unmapped {ref_file} {in_fastq}"
            do.run(cmd.format(**locals()), "Index STAR")
            if os.path.exists(out_dir):
                shutil.rmtree(out_dir)
            shutil.move(tx_out_dir, out_dir)
    return out_dir

def calculate_bismark_instances(cores, memory):
    """
    calculate number of parallel bismark instances to run, based on disussion here
    https://github.com/FelixKrueger/Bismark/issues/96
    cores and memory here are the maximum amounts available for us to use
    """
    BISMARK_CORES = 1
    BOWTIE_CORES_PER_INSTANCE = 2
    SAMTOOLS_CORES_PER_INSTANCE = 1
    CORES_PER_INSTANCE = BOWTIE_CORES_PER_INSTANCE + SAMTOOLS_CORES_PER_INSTANCE
    GENOME_MEMORY_GB = 12
    INSTANCE_MEMORY_GB = 10

    available_instance_memory = memory - GENOME_MEMORY_GB
    instances_in_memory = max(available_instance_memory / INSTANCE_MEMORY_GB, 1)

    available_instance_cores = cores - BISMARK_CORES
    instances_in_cores = max(available_instance_cores / CORES_PER_INSTANCE, 1)
    instances = int(min(instances_in_memory, instances_in_cores))
    logger.info(f"{cores} cores and {memory} memory are available. Spinning up {instances} instances of bismark.")
    return instances
