import os

from bcbio import utils
from bcbio.distributed import transaction
from bcbio.pipeline import datadict, config_utils
from bcbio.provenance import do
from bcbio import bam


def dedup_bismark(data):
    """Remove alignments to the same position in the genome from the Bismark
    mapping output using deduplicate_bismark
    """
    input_file = datadict.get_work_bam(data)
    input_file = bam.sort(input_file, datadict.get_config(data), order="queryname")
    sample_name = datadict.get_sample_name(data)
    output_dir = os.path.join(datadict.get_work_dir(data), 'dedup',
                              sample_name)
    output_dir = utils.safe_makedir(output_dir)

    input_file_name, input_file_extension = os.path.splitext(os.path.basename(
        input_file
    ))
    output_file = os.path.join(
        output_dir, f'{input_file_name}.deduplicated{input_file_extension}'
    )

    if utils.file_exists(output_file):
        data = datadict.set_work_bam(data, output_file)
        return [[data]]

    deduplicate_bismark = config_utils.get_program('deduplicate_bismark',
                                                   data['config'])
    command = f'{deduplicate_bismark} --output_dir {output_dir} {input_file}'
    with transaction.file_transaction(output_dir):
        do.run(command, 'remove deduplicate alignments')

    data = datadict.set_work_bam(data, output_file)
    return [[data]]
