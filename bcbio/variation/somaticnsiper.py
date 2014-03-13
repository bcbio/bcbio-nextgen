"""Provide support for SomaticSniper paired analysis."""

from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.variation.realign import has_aligned_reads
from bcbio.pipeline import config_utils
from bcbio.variation import bamprep
from bcbio.variation.vcfutils import write_empty_vcf, get_paired_bams


def somaticsniper_caller(align_bams, items, ref_file, assoc_files, region=None,
                       out_file=None):

    if out_file is None:
        out_file = "%s-paired-variants.vcf" % os.path.splitext(
            align_bams[0])[0]

    base_config = items[0]["config"]
    somaticsniper = config_utils.get_program("somaticsniper", base_config)

    tumor_bam, tumor_name, normal_bam, normal_name = get_paired_bams(
        align_bams, items)

    if not file_exists(out_file):
        #FIXME: A version check would be needed, but somaticsniper had
        # *no* official release with the fix we need

        with file_transaction(out_file) as tx_out:

            params = ("{somaticsniper} -f {ref_file} -n {normal_name} "
                      "-t {tumor_name} -F vcf {tumor_bam} {normal_bam} "
                      "{tx_out}")
            cmd = params.format(**locals())
            do.run(cmd, None, [do.file_exists(tx_out)])

    return out_file
