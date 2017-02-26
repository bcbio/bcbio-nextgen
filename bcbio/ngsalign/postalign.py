"""Perform streaming post-alignment preparation -- de-duplication and sorting.

Centralizes a pipelined approach to generating sorted, de-duplicated BAM output
from sequencer results.

samblaster: http://arxiv.org/pdf/1403.7486v1.pdf
biobambam bammarkduplicates: http://arxiv.org/abs/1306.0836
"""
import contextlib
import os

import toolz as tz

from bcbio import bam, broad, utils
from bcbio.bam import ref
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.log import logger
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do

@contextlib.contextmanager
def tobam_cl(data, out_file, is_paired=False):
    """Prepare command line for producing de-duplicated sorted output.

    - If no deduplication, sort and prepare a BAM file.
    - If paired, then use samblaster and prepare discordant outputs.
    - If unpaired, use biobambam's bammarkduplicates
    """
    do_dedup = _check_dedup(data)
    umi_consensus = dd.get_umi_consensus(data)
    with file_transaction(data, out_file) as tx_out_file:
        if not do_dedup:
            yield (sam_to_sortbam_cl(data, tx_out_file), tx_out_file)
        elif umi_consensus:
            yield (_sam_to_grouped_umi_cl(data, umi_consensus, tx_out_file), tx_out_file)
        elif is_paired and _need_sr_disc_reads(data) and not _too_many_contigs(dd.get_ref_file(data)):
            sr_file = "%s-sr.bam" % os.path.splitext(out_file)[0]
            disc_file = "%s-disc.bam" % os.path.splitext(out_file)[0]
            with file_transaction(data, sr_file) as tx_sr_file:
                with file_transaction(data, disc_file) as tx_disc_file:
                    yield (samblaster_dedup_sort(data, tx_out_file, tx_sr_file, tx_disc_file),
                           tx_out_file)
        else:
            yield (_biobambam_dedup_sort(data, tx_out_file), tx_out_file)

def _too_many_contigs(ref_file):
    """Check for more contigs than the maximum samblaster deduplication supports.
    """
    max_contigs = 32768
    return len(list(ref.file_contigs(ref_file))) >= max_contigs

def _need_sr_disc_reads(data):
    """Check if we need split and discordant reads in downstream processing.

    We use samblaster when needed and otherwise use an approach that does not
    extract these reads to be less resource intensive.
    """
    from bcbio import structural
    return "lumpy" in structural.get_svcallers(data)

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
    tmp_prefix = "%s-sorttmp" % utils.splitext_plus(tx_out_file)[0]
    tobam_cmd = ("{samtools} sort {sort_opt} -@ {cores} -m {mem} -T {tmp_prefix}-{dext} -o {out_file} -")
    # full BAM -- associate more memory and cores
    cores, mem = _get_cores_memory(data, downscale=2)
    sort_opt = "-n" if data.get("align_split") else ""
    dedup_cmd = tobam_cmd.format(out_file=tx_out_file, dext="full", **locals())
    # split and discordant BAMs -- give less memory/cores since smaller files
    sort_opt = ""
    cores, mem = _get_cores_memory(data, downscale=4)
    splitter_cmd = tobam_cmd.format(out_file=tx_sr_file, dext="spl", **locals())
    discordant_cmd = tobam_cmd.format(out_file=tx_disc_file, dext="disc", **locals())
    # samblaster 0.1.22 and better require the -M flag for compatibility with bwa-mem
    cmd = ("{samblaster} --addMateTags -M --splitterFile >({splitter_cmd}) --discordantFile >({discordant_cmd}) "
           "| {dedup_cmd}")
    return cmd.format(**locals())

def _biobambam_dedup_sort(data, tx_out_file):
    """Perform streaming deduplication and sorting with biobambam's bamsormadup
    """
    samtools = config_utils.get_program("samtools", data["config"])
    cores, mem = _get_cores_memory(data, downscale=2)
    tmp_file = "%s-sorttmp" % utils.splitext_plus(tx_out_file)[0]
    if data.get("align_split"):
        cmd = "{samtools} sort -n -@ {cores} -m {mem} -O bam -T {tmp_file}-namesort -o {tx_out_file} -"
    else:
        cmd = ("bamsormadup inputformat=sam threads={cores} tmpfile={tmp_file}-markdup "
               "SO=coordinate indexfilename={tx_out_file}.bai > {tx_out_file}")
    return cmd.format(**locals())

def _sam_to_grouped_umi_cl(data, umi_consensus, tx_out_file):
    """Mark duplicates on aligner output and convert to grouped UMIs by position.

    Works with either a separate umi_file or UMI embedded in the read names.
    """
    tmp_file = "%s-sorttmp" % utils.splitext_plus(tx_out_file)[0]
    jvm_opts = _get_fgbio_jvm_opts(data, os.path.dirname(tmp_file), 1)
    cores, mem = _get_cores_memory(data)
    cmd = ("bamsormadup tmpfile={tmp_file}-markdup inputformat=sam threads={cores} outputformat=bam "
           "level=0 SO=coordinate | ")
    # UMIs in a separate file
    if os.path.exists(umi_consensus):
        cmd += "fgbio {jvm_opts} AnnotateBamWithUmis -i /dev/stdin -f {umi_consensus} -o {tx_out_file}"
    # UMIs embedded in read name
    else:
        cmd += "umis bamtag - | samtools view -b > {tx_out_file}"
    return cmd.format(**locals())

def _get_fgbio_jvm_opts(data, tmpdir, scale_factor=None):
    cores, mem = _get_cores_memory(data)
    resources = config_utils.get_resources("fgbio", data["config"])
    jvm_opts = resources.get("jvm_opts", ["-Xms750m", "-Xmx4g"])
    if scale_factor and cores > scale_factor:
        jvm_opts = config_utils.adjust_opts(jvm_opts, {"algorithm": {"memory_adjust":
                                                                     {"direction": "increase",
                                                                      "magnitude": cores // scale_factor}}})
    jvm_opts += broad.get_default_jvm_opts(tmpdir)
    jvm_opts = " ".join(jvm_opts)
    return jvm_opts

def umi_consensus(data):
    """Convert UMI grouped reads into fastq pair for re-alignment.
    """
    align_bam = dd.get_work_bam(data)
    f1_out = "%s-cumi-1.fq.gz" % utils.splitext_plus(align_bam)[0]
    f2_out = "%s-cumi-2.fq.gz" % utils.splitext_plus(align_bam)[0]
    if not utils.file_uptodate(f1_out, align_bam):
        with file_transaction(data, f1_out, f2_out) as (tx_f1_out, tx_f2_out):
            jvm_opts = _get_fgbio_jvm_opts(data, os.path.dirname(tx_f1_out), 2)
            group_opts, cons_opts = _get_fgbio_options(data)
            cmd = ("unset JAVA_HOME && "
                   "fgbio {jvm_opts} GroupReadsByUmi {group_opts} -s adjacency -i {align_bam} | "
                   "fgbio {jvm_opts} CallMolecularConsensusReads {cons_opts} "
                   "-S queryname -i /dev/stdin -o /dev/stdout | "
                   "bamtofastq F={tx_f1_out} F2={tx_f2_out} gz=1")
            do.run(cmd.format(**locals()), "UMI consensus fastq generation")
    return f1_out, f2_out

def _get_fgbio_options(data):
    """Get adjustable, through resources, or default options for fgbio.
    """
    group_opts = ["--edits", "--min-map-q"]
    cons_opts = ["--min-reads", "--min-consensus-base-quality", "--min-input-base-quality"]
    defaults = {"--min-reads": "1",
                "--min-map-q": "1",
                "--min-consensus-base-quality": "13",
                "--min-input-base-quality": "2",
                "--edits": "1"}
    ropts = config_utils.get_resources("fgbio", data["config"]).get("options", [])
    assert len(ropts) % 2 == 0, "Expect even number of options for fgbio" % ropts
    defaults.update(dict(tz.partition(2, ropts)))
    group_out = " ".join(["%s %s" % (x, defaults[x]) for x in group_opts])
    cons_out = " ".join(["%s %s" % (x, defaults[x]) for x in cons_opts])
    return group_out, cons_out

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
