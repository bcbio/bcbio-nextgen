"""Task definitions for the Celery message queue (http://celeryproject.org/).
"""
import time

from celery.task import task

from bcbio.pipeline import sample, lane, qcsummary, toplevel, storage, shared, variation, validate
from bcbio import structural
from bcbio.variation import realign, genotype, ensemble, population, multi, recalibrate, vcfutils
from bcbio import chipseq

# Global configuration for tasks in the main celeryconfig module
import celeryconfig

@task(ignore_results=True, queue="toplevel")
def analyze_and_upload(*args):
    """Run full analysis and upload results to Galaxy instance.

    Workers need to run on the machine with Galaxy installed for upload,
    but the actual processing can be distributed to multiple nodes.
    """
    config_file = celeryconfig.BCBIO_CONFIG_FILE
    remote_info = args[0]
    toplevel.analyze_and_upload(remote_info, config_file)

@task(ignore_results=True, queue="storage")
def long_term_storage(*args):
    config_file = celeryconfig.BCBIO_CONFIG_FILE
    remote_info = args[0]
    storage.long_term_storage(remote_info, config_file)

@task
def process_lane(*args):
    return lane.process_lane(*args)

@task
def trim_lane(*args):
    return lane.trim_lane(*args)

@task
def process_alignment(*args):
    return lane.process_alignment(*args)

@task
def merge_sample(*args):
    return sample.merge_sample(*args)

@tasks
def delayed_bam_merge(*args):
    return sample.delayed_bam_merge(*args)

@task
def prep_recal(*args):
    return recalibrate.prep_recal(*args)

@task
def write_recal_bam(*args):
    return recalibrate.write_recal_bam(*args)

@task
def realign_sample(*args):
    return realign.realign_sample(*args)

@task
def pipeline_summary(*args):
    return qcsummary.pipeline_summary(*args)

@task
def generate_transcript_counts(*args):
    return sample.generate_transcript_counts(*args)

@task
def split_variants_by_sample(*args):
    return multi.split_variants_by_sample(*args)

@task
def postprocess_variants(*args):
    return variation.postprocess_variants(*args)

@task
def generate_bigwig(*args):
    return sample.generate_bigwig(*args)

@task
def combine_bam(*args):
    return shared.combine_bam(*args)

@task
def variantcall_sample(*args):
    return genotype.variantcall_sample(*args)

@task
def combine_variant_files(*args):
    return vcfutils.combine_variant_files(*args)

@task
def concat_variant_files(*args):
    return vcfutils.concat_variant_files(*args)

@task
def detect_sv(*args):
    return structural.detect_sv(*args)

@task
def combine_calls(*args):
    return ensemble.combine_calls(*args)

@task
def prep_gemini_db(*args):
    return population.prep_gemini_db(*args)

@task
def compare_to_rm(*args):
    return validate.compare_to_rm(*args)

@task
def test(x):
    print x
    time.sleep(5)
    return x

@task
def clean_chipseq_aligment(*args):
    return chipseq.clean_chipseq_alignment(*args)
