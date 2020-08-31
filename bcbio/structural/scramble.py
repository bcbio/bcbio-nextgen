"""Mobile element insertion detection using SCRAMble (Soft Clipped Read Alignment Mapper)

https://github.com/GeneDx/scramble
https://doi.org/10.1038/s41436-020-0749-x
"""
import os

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do


def run(items):
    for index, data in enumerate(items):
        work_dir = os.path.join(dd.get_work_dir(data), 'structural', dd.get_sample_name(data))
        output_clusters_file_path = os.path.join(work_dir, 'clusters.txt')

        if not utils.file_exists(output_clusters_file_path):
            with file_transaction(output_clusters_file_path) as temp_output_file_path:
                do.run(f'cluster_identifier {dd.get_work_bam(data)} > {temp_output_file_path}')

        items[index].get('sv', []).append({'variantcaller': 'scramble',
                                           'clusters_file': output_clusters_file_path})
    return items
