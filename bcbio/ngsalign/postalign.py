"""Perform streaming post-alignment preparation -- de-duplication and sorting.

Centralizes a pipelined approach to generating sorted, de-duplicated BAM output
from sequencer results.

sambamba: https://github.com/lomereiter/sambamba
samblaster: http://arxiv.org/pdf/1403.7486v1.pdf
biobambam bammarkduplicates: http://arxiv.org/abs/1306.0836
"""
import contextlib
import os

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import config_utils

@contextlib.contextmanager
def tobam_cl(data, out_file, is_paired=False):
    """Prepare command line for producing de-duplicated sorted output.

    - If no deduplication, sort and prepare a BAM file.
    - If paired, then use samblaster and prepare discordant outputs.
    - If unpaired, use biobambam's bammarkduplicates
    """
    do_dedup = _check_dedup(data)
    with utils.curdir_tmpdir(data) as tmpdir:
        with file_transaction(out_file) as tx_out_file:
            if not do_dedup:
                yield (_sam_to_sortbam_cl(data, tmpdir, tx_out_file), tx_out_file)
            elif is_paired:
                sr_file = "%s-sr.bam" % os.path.splitext(out_file)[0]
                disc_file = "%s-disc.bam" % os.path.splitext(out_file)[0]
                with file_transaction(sr_file) as tx_sr_file:
                    with file_transaction(disc_file) as tx_disc_file:
                        yield (samblaster_dedup_sort(data, tmpdir, tx_out_file, tx_sr_file, tx_disc_file),
                               tx_out_file)
            else:
                yield (_biobambam_dedup_sort(data, tmpdir, tx_out_file), tx_out_file)

def _get_cores_memory(data, downscale=2):
    """Retrieve cores and memory, using samtools as baseline.

    For memory, scaling down because we share with alignment and de-duplication.
    """
    resources = config_utils.get_resources("samtools", data["config"])
    num_cores = data["config"]["algorithm"].get("num_cores", 1)
    max_mem = config_utils.adjust_memory(resources.get("memory", "2G"),
                                         downscale, "decrease").upper()
    return num_cores, max_mem

def _sam_to_sortbam_cl(data, tmpdir, tx_out_file):
    """Convert to sorted BAM output with sambamba.
    """
    samtools = config_utils.get_program("samtools", data["config"])
    sambamba = config_utils.get_program("sambamba", data["config"])
    cores, mem = _get_cores_memory(data, downscale=3)
    return ("{samtools} view -b -S -u - | "
            "{sambamba} sort -t {cores} -m {mem} "
            "--tmpdir {tmpdir} -o {tx_out_file} /dev/stdin".format(**locals()))

def samblaster_dedup_sort(data, tmpdir, tx_out_file, tx_sr_file, tx_disc_file):
    """Deduplicate and sort with samblaster, produces split read and discordant pair files.
    """
    sambamba = config_utils.get_program("sambamba", data["config"])
    samblaster = config_utils.get_program("samblaster", data["config"])
    samtools = config_utils.get_program("samtools", data["config"])
    cores, mem = _get_cores_memory(data, downscale=3)
    for dname in ["spl", "disc", "full"]:
        utils.safe_makedir(os.path.join(tmpdir, dname))
    tobam_cmd = ("{samtools} view -S -u /dev/stdin | "
                 "{sambamba} sort -t {cores} -m {mem} --tmpdir {tmpdir}/{dext} "
                 "-o {out_file} /dev/stdin")
    splitter_cmd = tobam_cmd.format(out_file=tx_sr_file, dext="spl", **locals())
    discordant_cmd = tobam_cmd.format(out_file=tx_disc_file, dext="disc", **locals())
    dedup_cmd = tobam_cmd.format(out_file=tx_out_file, dext="full", **locals())
    cmd = ("{samblaster} --splitterFile >({splitter_cmd}) --discordantFile >({discordant_cmd}) "
           "| {dedup_cmd}")
    return cmd.format(**locals())

def _biobambam_dedup_sort(data, tmpdir, tx_out_file):
    """Perform streaming deduplication and sorting with biobambam's bammarkduplicates2.
    """
    samtools = config_utils.get_program("samtools", data["config"])
    sambamba = config_utils.get_program("sambamba", data["config"])
    bammarkduplicates = config_utils.get_program("bammarkduplicates", data["config"])
    base_tmp = os.path.join(tmpdir, os.path.splitext(os.path.basename(tx_out_file))[0])
    cores, mem = _get_cores_memory(data, downscale=3)
    sort2_tmpdir = utils.safe_makedir(os.path.join(tmpdir, "sort2"))
    return ("{samtools} view -b -S -u - |"
            "{samtools} sort -n -o -@ {cores} -m {mem} - {base_tmp}-sort | "
            "{bammarkduplicates} tmpfile={base_tmp}-markdup "
            "markthreads={cores} level=0 | "
            "{sambamba} sort -t {cores} -m {mem} --tmpdir {sort2_tmpdir} "
            "-o {tx_out_file} /dev/stdin").format(**locals())

def _check_dedup(data):
    """Check configuration for de-duplication, handling back compatibility.
    """
    dup_param = utils.get_in(data, ("config", "algorithm", "mark_duplicates"), True)
    if dup_param and isinstance(dup_param, basestring):
        logger.info("Warning: bcbio no longer support explicit setting of mark_duplicate algorithm. "
                    "Using best-practice choice based on input data.")
        dup_param = True
    return dup_param
