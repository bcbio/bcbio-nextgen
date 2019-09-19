"""Perform streaming post-alignment preparation -- de-duplication and sorting.

Centralizes a pipelined approach to generating sorted, de-duplicated BAM output
from sequencer results.

samblaster: http://arxiv.org/pdf/1403.7486v1.pdf
biobambam bammarkduplicates: http://arxiv.org/abs/1306.0836
"""
import contextlib
import math
import os

import toolz as tz

from bcbio import bam, broad, utils
from bcbio.bam import ref
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.log import logger
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import coverage

import six


pysam = utils.LazyImport("pysam")

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
    tobam_cmd = ("{samtools} sort {sort_opt} -@ {cores} -m {mem} -T {tmp_prefix}-{dext} {out_file} -")
    # full BAM -- associate more memory and cores
    cores, mem = _get_cores_memory(data, downscale=2)
    # Potentially downsample to maximum coverage here if not splitting and whole genome sample
    ds_cmd = None if data.get("align_split") else bam.get_maxcov_downsample_cl(data, "samtools")
    sort_opt = "-n" if data.get("align_split") and dd.get_mark_duplicates(data) else ""
    if ds_cmd:
        dedup_cmd = "%s %s > %s" % (tobam_cmd.format(out_file="", dext="full", **locals()), ds_cmd, tx_out_file)
    else:
        dedup_cmd = tobam_cmd.format(out_file="-o %s" % tx_out_file, dext="full", **locals())
    # split and discordant BAMs -- give less memory/cores since smaller files
    sort_opt = ""
    cores, mem = _get_cores_memory(data, downscale=4)
    splitter_cmd = tobam_cmd.format(out_file="-o %s" % tx_sr_file, dext="spl", **locals())
    discordant_cmd = tobam_cmd.format(out_file="-o %s" % tx_disc_file, dext="disc", **locals())
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
        sort_opt = "-n" if data.get("align_split") and _check_dedup(data) else ""
        cmd = "{samtools} sort %s -@ {cores} -m {mem} -O bam -T {tmp_file}-namesort -o {tx_out_file} -" % sort_opt
    else:
        # scale core usage to avoid memory issues with larger WGS samples
        cores = max(1, int(math.ceil(cores * 0.75)))
        ds_cmd = bam.get_maxcov_downsample_cl(data, "bamsormadup")
        bamsormadup = config_utils.get_program("bamsormadup", data)
        cmd = ("{bamsormadup} inputformat=sam threads={cores} tmpfile={tmp_file}-markdup "
               "SO=coordinate %s > {tx_out_file}" % ds_cmd)
    return cmd.format(**locals())

def _sam_to_grouped_umi_cl(data, umi_consensus, tx_out_file):
    """Mark duplicates on aligner output and convert to grouped UMIs by position.

    Works with either a separate umi_file or UMI embedded in the read names.
    """
    tmp_file = "%s-sorttmp" % utils.splitext_plus(tx_out_file)[0]
    jvm_opts = _get_fgbio_jvm_opts(data, os.path.dirname(tmp_file), 1)
    cores, mem = _get_cores_memory(data)
    bamsormadup = config_utils.get_program("bamsormadup", data)
    cmd = ("{bamsormadup} tmpfile={tmp_file}-markdup inputformat=sam threads={cores} outputformat=bam "
           "level=0 SO=coordinate | ")
    # UMIs in a separate file
    if os.path.exists(umi_consensus) and os.path.isfile(umi_consensus):
        cmd += "fgbio {jvm_opts} AnnotateBamWithUmis -i /dev/stdin -f {umi_consensus} -o {tx_out_file}"
    # UMIs embedded in read name
    else:
        cmd += ("%s %s bamtag - | samtools view -b > {tx_out_file}" %
                (utils.get_program_python("umis"),
                 config_utils.get_program("umis", data["config"])))
    return cmd.format(**locals())

def _get_fgbio_jvm_opts(data, tmpdir, scale_factor=None):
    cores, mem = _get_cores_memory(data)
    resources = config_utils.get_resources("fgbio", data["config"])
    jvm_opts = resources.get("jvm_opts", ["-Xms750m", "-Xmx4g"])
    if scale_factor and cores > scale_factor:
        jvm_opts = config_utils.adjust_opts(jvm_opts, {"algorithm": {"memory_adjust":
                                                                     {"direction": "increase",
                                                                      "magnitude": cores // scale_factor}}})
    jvm_opts += broad.get_default_jvm_opts()
    jvm_opts = " ".join(jvm_opts)
    return jvm_opts + " --tmp-dir %s" % tmpdir

def _estimate_fgbio_defaults(avg_coverage):
    """Provide fgbio defaults based on input sequence depth and coverage.

    For higher depth/duplication we want to use `--min-reads` to allow
    consensus calling in the duplicates:

    https://fulcrumgenomics.github.io/fgbio/tools/latest/CallMolecularConsensusReads.html

    If duplicated adjusted depth leaves a coverage of 800x or higher
    (giving us ~4 reads at 0.5% detection frequency),
    then we use `--min-reads 2`, otherwise `--min-reads 1`
    """
    out = {}
    if avg_coverage >= 800:
        out["--min-reads"] = 2
    else:
        out["--min-reads"] = 1
    return out

def correct_umis(data):
    """Correct umis against the whitelist in correct_umi_file

    http://fulcrumgenomics.github.io/fgbio/tools/latest/CorrectUmis.html
    """
    input_bam = dd.get_work_bam(data)
    output_bam = os.path.join(utils.safe_makedir(os.path.join(os.getcwd(),
                              "align", dd.get_sample_name(data))),
                              "%s-umis_corrected%s" % utils.splitext_plus(os.path.basename(input_bam)))
    jvm_opts = _get_fgbio_jvm_opts(data, os.path.dirname(output_bam), 2)
    # Improve speeds by avoiding compression read/write bottlenecks
    io_opts = "--async-io=true --compression=0"
    umis_whitelist = tz.get_in(["config", "algorithm", "correct_umis"], data)
    umi_method, umi_tag = _check_umi_type(input_bam)

    cmd = ("unset JAVA_HOME && "
           "fgbio {jvm_opts} {io_opts} CorrectUmis "
           "-t {umi_tag} -m 3 -d 1 -x "
           "-U {umis_whitelist} "
           "-i {input_bam} -o {output_bam}")
    do.run(cmd.format(**locals()), "Correcting UMIs")
    return output_bam

def umi_consensus(data):
    """Convert UMI grouped reads into fastq pair for re-alignment.
    """
    align_bam = dd.get_work_bam(data)
    umi_method, umi_tag = _check_umi_type(align_bam)
    f1_out = "%s-cumi-1.fq.gz" % utils.splitext_plus(align_bam)[0]
    f2_out = "%s-cumi-2.fq.gz" % utils.splitext_plus(align_bam)[0]
    avg_coverage = coverage.get_average_coverage("rawumi", dd.get_variant_regions(data), data)
    if not utils.file_uptodate(f1_out, align_bam):
        with file_transaction(data, f1_out, f2_out) as (tx_f1_out, tx_f2_out):
            jvm_opts = _get_fgbio_jvm_opts(data, os.path.dirname(tx_f1_out), 2)
            # Improve speeds by avoiding compression read/write bottlenecks
            io_opts = "--async-io=true --compression=0"
            est_options = _estimate_fgbio_defaults(avg_coverage)
            group_opts, cons_opts, filter_opts = _get_fgbio_options(data, est_options, umi_method)
            cons_method = "CallDuplexConsensusReads" if umi_method == "paired" else "CallMolecularConsensusReads"
            tempfile = "%s-bamtofastq-tmp" % utils.splitext_plus(f1_out)[0]
            ref_file = dd.get_ref_file(data)
            cmd = ("unset JAVA_HOME && "
                   "fgbio {jvm_opts} {io_opts} GroupReadsByUmi {group_opts} -t {umi_tag} -s {umi_method} "
                   "-i {align_bam} | "
                   "fgbio {jvm_opts} {io_opts} {cons_method} {cons_opts} --sort-order=:none: "
                   "-i /dev/stdin -o /dev/stdout | "
                   "fgbio {jvm_opts} {io_opts} FilterConsensusReads {filter_opts} -r {ref_file} "
                   "-i /dev/stdin -o /dev/stdout | "
                   "bamtofastq collate=1 T={tempfile} F={tx_f1_out} F2={tx_f2_out} tags=cD,cM,cE gz=1")
            do.run(cmd.format(**locals()), "UMI consensus fastq generation")
    return f1_out, f2_out, avg_coverage

def _check_umi_type(bam_file):
    """Determine the type of UMI from BAM tags: standard or paired.
    """
    with pysam.Samfile(bam_file, "rb") as in_bam:
        for read in in_bam:
            cur_umi = None
            for tag in ["RX", "XC"]:
                try:
                    cur_umi = read.get_tag(tag)
                    break
                except KeyError:
                    pass
            if cur_umi:
                if "-" in cur_umi and len(cur_umi.split("-")) == 2:
                    return "paired", tag
                else:
                    return "adjacency", tag

def _get_fgbio_options(data, estimated_defaults, umi_method):
    """Get adjustable, through resources, or default options for fgbio.
    """
    group_opts = ["--edits", "--min-map-q"]
    cons_opts = ["--min-input-base-quality"]
    if umi_method != "paired":
        cons_opts += ["--min-reads", "--max-reads"]
    filter_opts = ["--min-reads", "--min-base-quality", "--max-base-error-rate"]
    defaults = {"--min-reads": "1",
                "--max-reads": "100000",
                "--min-map-q": "1",
                "--min-base-quality": "13",
                "--max-base-error-rate": "0.1",
                "--min-input-base-quality": "2",
                "--edits": "1"}
    defaults.update(estimated_defaults)
    ropts = config_utils.get_resources("fgbio", data["config"]).get("options", [])
    assert len(ropts) % 2 == 0, "Expect even number of options for fgbio" % ropts
    ropts = dict(tz.partition(2, ropts))
    # Back compatibility for older base quality settings
    if "--min-consensus-base-quality" in ropts:
        ropts["--min-base-quality"] = ropts.pop("--min-consensus-base-quality")
    defaults.update(ropts)
    group_out = " ".join(["%s=%s" % (x, defaults[x]) for x in group_opts])
    cons_out = " ".join(["%s=%s" % (x, defaults[x]) for x in cons_opts])
    filter_out = " ".join(["%s=%s" % (x, defaults[x]) for x in filter_opts])
    if umi_method != "paired":
        cons_out += " --output-per-base-tags=false"
    return group_out, cons_out, filter_out

def _check_dedup(data):
    """Check configuration for de-duplication.

    Defaults to no de-duplication for RNA-seq and small RNA, the
    back compatible default. Allow overwriting with explicit
    `mark_duplicates: true` setting.
    Also defaults to false for no alignment inputs.
    """
    if dd.get_analysis(data).lower() in ["rna-seq", "smallrna-seq"] or not dd.get_aligner(data):
        dup_param = utils.get_in(data, ("config", "algorithm", "mark_duplicates"), False)
    else:
        dup_param = utils.get_in(data, ("config", "algorithm", "mark_duplicates"), True)
    if dup_param and isinstance(dup_param, six.string_types):
        logger.info("Warning: bcbio no longer support explicit setting of mark_duplicate algorithm. "
                    "Using best-practice choice based on input data.")
        dup_param = True
    return dup_param

def dedup_bam(in_bam, data):
    """Perform non-stream based deduplication of BAM input files using biobambam.
    """
    if _check_dedup(data):
        out_file = os.path.join(utils.safe_makedir(os.path.join(os.getcwd(), "align", dd.get_sample_name(data))),
                                "%s-dedup%s" % utils.splitext_plus(os.path.basename(in_bam)))
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
