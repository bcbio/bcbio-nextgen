"""Perform streaming post-alignment preparation -- de-duplication and sorting.

Centralizes a pipelined approach to generating sorted, de-duplicated BAM output
from sequencer results.

sambamba: https://github.com/lomereiter/sambamba
samblaster: http://arxiv.org/pdf/1403.7486v1.pdf
biobambam bammarkduplicates: http://arxiv.org/abs/1306.0836
"""
import contextlib
from distutils.version import LooseVersion
import os

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.log import logger
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do, programs

@contextlib.contextmanager
def tobam_cl(data, out_file, is_paired=False):
    """Prepare command line for producing de-duplicated sorted output.

    - If no deduplication, sort and prepare a BAM file.
    - If paired, then use samblaster and prepare discordant outputs.
    - If unpaired, use biobambam's bammarkduplicates
    """
    do_dedup = _check_dedup(data)
    with file_transaction(data, out_file) as tx_out_file:
        if not do_dedup:
            yield (sam_to_sortbam_cl(data, tx_out_file), tx_out_file)
        elif is_paired:
            sr_file = "%s-sr.bam" % os.path.splitext(out_file)[0]
            disc_file = "%s-disc.bam" % os.path.splitext(out_file)[0]
            with file_transaction(data, sr_file) as tx_sr_file:
                with file_transaction(data, disc_file) as tx_disc_file:
                    yield (samblaster_dedup_sort(data, tx_out_file, tx_sr_file, tx_disc_file),
                           tx_out_file)
        else:
            yield (_biobambam_dedup_sort(data, tx_out_file), tx_out_file)

def _get_cores_memory(data, downscale=2):
    """Retrieve cores and memory, using samtools as baseline.

    For memory, scaling down because we share with alignment and de-duplication.
    """
    resources = config_utils.get_resources("samtools", data["config"])
    num_cores = data["config"]["algorithm"].get("num_cores", 1)
    max_mem = config_utils.adjust_memory(resources.get("memory", "2G"),
                                         downscale, "decrease").upper()
    return num_cores, max_mem

def sam_to_sortbam_cl(data, tx_out_file, name_sort=False):
    """Convert to sorted BAM output.

    Set name_sort to True to sort reads by queryname
    """
    samtools = config_utils.get_program("samtools", data["config"])
    cores, mem = _get_cores_memory(data, downscale=2)
    tmp_file = "%s-sorttmp" % utils.splitext_plus(tx_out_file)[0]
    sort_flag = "-n" if name_sort else ""
    return ("{samtools} sort -@ {cores} -m {mem} {sort_flag} "
            "-T {tmp_file} -o {tx_out_file} /dev/stdin".format(**locals()))

def samblaster_dedup_sort(data, tx_out_file, tx_sr_file, tx_disc_file):
    """Deduplicate and sort with samblaster, produces split read and discordant pair files.
    """
    samblaster = config_utils.get_program("samblaster", data["config"])
    samtools = config_utils.get_program("samtools", data["config"])
    cores, mem = _get_cores_memory(data, downscale=3)
    tmp_prefix = "%s-sorttmp" % utils.splitext_plus(tx_out_file)[0]
    for ext in ["spl", "disc", "full"]:
        utils.safe_makedir("%s-%s" % (tmp_prefix, ext))
    sort_opt = "-N" if data.get("align_split") else ""
    full_tobam_cmd = ("samtools view -b -S -u - | "
                      "sambamba sort {sort_opt} -t {cores} -m {mem} "
                      "--tmpdir {tmp_prefix}-{dext} -o {out_file} /dev/stdin")
    tobam_cmd = ("{samtools} sort -@ {cores} -m {mem} "
                 "-T {tmp_prefix}-{dext} -o {out_file} /dev/stdin")
    # samblaster 0.1.22 and better require the -M flag for compatibility with bwa-mem
    # https://github.com/GregoryFaust/samblaster/releases/tag/v.0.1.22
    if LooseVersion(programs.get_version_manifest("samblaster", data=data, required=True)) >= LooseVersion("0.1.22"):
        opts = "-M"
    else:
        opts = ""
    splitter_cmd = tobam_cmd.format(out_file=tx_sr_file, dext="spl", **locals())
    discordant_cmd = tobam_cmd.format(out_file=tx_disc_file, dext="disc", **locals())
    dedup_cmd = full_tobam_cmd.format(out_file=tx_out_file, dext="full", **locals())
    cmd = ("{samblaster} --addMateTags {opts} --splitterFile >({splitter_cmd}) --discordantFile >({discordant_cmd}) "
           "| {dedup_cmd}")
    return cmd.format(**locals())

def _biobambam_dedup_sort(data, tx_out_file):
    """Perform streaming deduplication and sorting with biobambam's bammarkduplicates2.
    """
    samtools = config_utils.get_program("samtools", data["config"])
    bammarkduplicates = config_utils.get_program("bammarkduplicates", data["config"])
    cores, mem = _get_cores_memory(data, downscale=2)
    tmp_file = "%s-sorttmp" % utils.splitext_plus(tx_out_file)[0]
    if data.get("align_split"):
        cmd = "{samtools} sort -n -@ {cores} -m {mem} -O bam -T {tmp_file}-namesort -o {tx_out_file} -"
    else:
        cmd = ("{samtools} sort -n -@ {cores} -m {mem} -O bam -T {tmp_file}-namesort - | "
               "{bammarkduplicates} tmpfile={tmp_file}-markdup "
               "markthreads={cores} level=0 | "
               "{samtools} sort -@ {cores} -m {mem} -T {tmp_file}-finalsort "
               "-o {tx_out_file} /dev/stdin")
    return cmd.format(**locals())

def _check_dedup(data):
    """Check configuration for de-duplication, handling back compatibility.
    """
    dup_param = utils.get_in(data, ("config", "algorithm", "mark_duplicates"), True)
    if dup_param and isinstance(dup_param, basestring):
        logger.info("Warning: bcbio no longer support explicit setting of mark_duplicate algorithm. "
                    "Using best-practice choice based on input data.")
        dup_param = True
    return dup_param

def dedup_bam(in_bam, data):
    """Perform non-stream based deduplication of BAM input files using biobambam.
    """
    if _check_dedup(data):
        out_file = "%s-dedup%s" % utils.splitext_plus(in_bam)
        if not utils.file_exists(out_file):
            with tx_tmpdir(data) as tmpdir:
                with file_transaction(data, out_file) as tx_out_file:
                    bammarkduplicates = config_utils.get_program("bammarkduplicates", data["config"])
                    base_tmp = os.path.join(tmpdir, os.path.splitext(os.path.basename(tx_out_file))[0])
                    cores, mem = _get_cores_memory(data, downscale=2)
                    cmd = ("{bammarkduplicates} tmpfile={base_tmp}-markdup "
                           "markthreads={cores} I={in_bam} O={tx_out_file}")
                    do.run(cmd.format(**locals()), "De-duplication with biobambam")
        bam.index(out_file, data["config"])
        return out_file
    else:
        return in_bam
