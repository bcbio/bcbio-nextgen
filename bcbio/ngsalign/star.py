import os
import sys
import shutil
import subprocess
import contextlib
from collections import namedtuple

import bcbio.bed as bed
from bcbio.pipeline import config_utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.utils import (safe_makedir, file_exists, is_gzipped)
from bcbio.provenance import do
from bcbio import utils

from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.ngsalign import postalign
from bcbio.bam import fastq

CLEANUP_FILES = ["Aligned.out.sam", "Log.out", "Log.progress.out"]
ALIGN_TAGS = ["NH", "HI", "NM", "MD", "AS"]


def align(fastq_file, pair_file, ref_file, names, align_dir, data):
    if not ref_file:
        logger.error("STAR index not found. We don't provide the STAR indexes "
                     "by default because they are very large. You can install "
                     "the index for your genome with: bcbio_nextgen.py upgrade "
                     "--aligners star --genomes genome-build-name --data")
        sys.exit(1)

    max_hits = 10
    srna = True if data["analysis"].lower().startswith("smallrna-seq") else False
    srna_opts = ""
    if srna:
        max_hits = 1000
        srna_opts = "--alignIntronMax 1"
    config = data["config"]
    star_dirs = _get_star_dirnames(align_dir, data, names)
    if file_exists(star_dirs.final_out):
        data = _update_data(star_dirs.final_out, star_dirs.out_dir, names, data)
        return data

    star_path = config_utils.get_program("STAR", config)
    def _unpack_fastq(f):
        """Use process substitution instead of readFilesCommand for gzipped inputs.

        Prevents issues on shared filesystems that don't support FIFO:
        https://github.com/alexdobin/STAR/issues/143
        """
        if f and is_gzipped(f):
            return "<(gunzip -c %s)" % f
        else:
            return f
    fastq_files = (" ".join([_unpack_fastq(fastq_file), _unpack_fastq(pair_file)])
                   if pair_file else _unpack_fastq(fastq_file))
    num_cores = dd.get_num_cores(data)
    gtf_file = dd.get_transcriptome_gtf(data)
    if not gtf_file:
        gtf_file = dd.get_gtf_file(data)
    if ref_file.endswith("chrLength"):
        ref_file = os.path.dirname(ref_file)

    with file_transaction(data, align_dir) as tx_align_dir:
        tx_star_dirnames = _get_star_dirnames(tx_align_dir, data, names)
        tx_out_dir, tx_out_file, tx_out_prefix, tx_final_out = tx_star_dirnames
        safe_makedir(tx_align_dir)
        safe_makedir(tx_out_dir)
        cmd = ("{star_path} --genomeDir {ref_file} --readFilesIn {fastq_files} "
               "--runThreadN {num_cores} --outFileNamePrefix {tx_out_prefix} "
               "--outReadsUnmapped Fastx --outFilterMultimapNmax {max_hits} "
               "--outStd BAM_Unsorted {srna_opts} "
               "--limitOutSJcollapsed 2000000 "
               "--outSAMtype BAM Unsorted "
               "--outSAMmapqUnique 60 "
               "--outSAMunmapped Within --outSAMattributes %s " % " ".join(ALIGN_TAGS))
        cmd += _add_sj_index_commands(fastq_file, ref_file, gtf_file) if not srna else ""
        cmd += _read_group_option(names)
        if dd.get_fusion_caller(data):
            if "arriba" in dd.get_fusion_caller(data):
                cmd += (
                    "--chimSegmentMin 10 --chimOutType WithinBAM SoftClip Junctions "
                    "--chimJunctionOverhangMin 10 --chimScoreMin 1 --chimScoreDropMax 30 "
                    "--chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 "
                    "--alignSJstitchMismatchNmax 5 -1 5 5 "
                    "--chimSegmentReadGapMax 3 ")
            else: 
                cmd += (" --chimSegmentMin 12 --chimJunctionOverhangMin 12 "
                    "--chimScoreDropMax 30 --chimSegmentReadGapMax 5 "
                    "--chimScoreSeparation 5 ")
                if "oncofuse" in dd.get_fusion_caller(data):
                    cmd += "--chimOutType Junctions "
                else:
                    cmd += "--chimOutType WithinBAM "
        strandedness = utils.get_in(data, ("config", "algorithm", "strandedness"),
                                    "unstranded").lower()
        if strandedness == "unstranded" and not srna:
            cmd += " --outSAMstrandField intronMotif "
        if not srna:
            cmd += " --quantMode TranscriptomeSAM "

        resources = config_utils.get_resources("star", data["config"])
        if resources.get("options", []):
            cmd += " " + " ".join([str(x) for x in resources.get("options", [])])
        cmd += " | " + postalign.sam_to_sortbam_cl(data, tx_final_out)
        cmd += " > {tx_final_out} "
        run_message = "Running STAR aligner on %s and %s" % (fastq_file, ref_file)
        do.run(cmd.format(**locals()), run_message, None)

    data = _update_data(star_dirs.final_out, star_dirs.out_dir, names, data)
    return data


StarOutDirs = namedtuple(
    'StarOutDirs',
    ['out_dir', 'out_file', 'out_prefix', 'final_out']
)


def _get_star_dirnames(align_dir, data, names):
    ALIGNED_OUT_FILE = "Aligned.out.sam"
    out_prefix = os.path.join(align_dir, dd.get_lane(data))
    out_file = out_prefix + ALIGNED_OUT_FILE
    out_dir = os.path.join(align_dir, "%s_star" % dd.get_lane(data))
    final_out = os.path.join(out_dir, "{0}.bam".format(names["sample"]))
    return StarOutDirs(out_dir, out_file, out_prefix, final_out)


def _add_sj_index_commands(fq1, ref_file, gtf_file):
    """
    newer versions of STAR can generate splice junction databases on thephfly
    this is preferable since we can tailor it to the read lengths
    """
    if _has_sj_index(ref_file):
        return ""
    else:
        rlength = fastq.estimate_maximum_read_length(fq1)
        cmd = " --sjdbGTFfile %s " % gtf_file
        cmd += " --sjdbOverhang %s " % str(rlength - 1)
        return cmd

def _has_sj_index(ref_file):
    """this file won't exist if we can do on the fly splice junction indexing"""
    return (file_exists(os.path.join(ref_file, "sjdbInfo.txt")) and
            (file_exists(os.path.join(ref_file, "transcriptInfo.tab"))))

def _update_data(align_file, out_dir, names, data):
    data = dd.set_work_bam(data, align_file)
    data = dd.set_align_bam(data, align_file)
    transcriptome_file = _move_transcriptome_file(out_dir, names)
    data = dd.set_transcriptome_bam(data, transcriptome_file)
    sjfile = get_splicejunction_file(out_dir, data)
    if sjfile:
        data = dd.set_starjunction(data, sjfile)
        sjbed = junction2bed(sjfile)
        data = dd.set_junction_bed(data, sjbed)
    sjchimfile = get_chimericjunction_file(out_dir, data)
    data = dd.set_chimericjunction(data, sjchimfile)
    return data

def _move_transcriptome_file(out_dir, names):
    out_file = os.path.join(out_dir, "{0}.transcriptome.bam".format(names["sample"]))
    star_file = os.path.join(out_dir, os.pardir,
                            "{0}Aligned.toTranscriptome.out.bam".format(names["lane"]))
    # if the out_file or the star_file doesn't exist, we didn't run the
    # transcriptome mapping
    if not file_exists(out_file):
        if not file_exists(star_file):
            return None
        else:
            shutil.move(star_file, out_file)
    return out_file

def _read_group_option(names):
    rg_id = names["rg"]
    rg_sample = names["sample"]
    rg_library = names["pl"]
    rg_platform_unit = names["pu"]
    rg_lb = ("LB:%s " % names.get("lb")) if names.get("lb") else ""

    return (" --outSAMattrRGline ID:{rg_id} PL:{rg_library} "
            "PU:{rg_platform_unit} SM:{rg_sample} {rg_lb}").format(**locals())

def remap_index_fn(ref_file):
    """Map sequence references to equivalent star indexes
    """
    return os.path.join(os.path.dirname(os.path.dirname(ref_file)), "star")

def index(ref_file, out_dir, data):
    """Create a STAR index in the defined reference directory.
    """
    (ref_dir, local_file) = os.path.split(ref_file)
    gtf_file = dd.get_gtf_file(data)
    if not utils.file_exists(gtf_file):
        raise ValueError("%s not found, could not create a star index." % (gtf_file))
    if not utils.file_exists(out_dir):
        with tx_tmpdir(data, os.path.dirname(out_dir)) as tx_out_dir:
            num_cores = dd.get_cores(data)
            cmd = ("STAR --genomeDir {tx_out_dir} --genomeFastaFiles {ref_file} "
                   "--runThreadN {num_cores} "
                   "--runMode genomeGenerate --sjdbOverhang 99 --sjdbGTFfile {gtf_file}")
            do.run(cmd.format(**locals()), "Index STAR")
            if os.path.exists(out_dir):
                shutil.rmtree(out_dir)
            shutil.move(tx_out_dir, out_dir)
    return out_dir

def get_star_version(data):
    star_path = config_utils.get_program("STAR", dd.get_config(data))
    cmd = "%s --version" % star_path
    subp = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            shell=True)
    with contextlib.closing(subp.stdout) as stdout:
        for line in stdout:
            if "STAR_" in line:
                version = line.split("STAR_")[1].strip()
    return version

def get_chimericjunction_file(out_dir, data):
    """
    locate the chimeric splice junction file starting from the alignment directory
    """
    samplename = dd.get_sample_name(data)
    sjfile = os.path.join(out_dir, os.pardir, f"{samplename}Chimeric.out.junction")
    if file_exists(sjfile):
        return sjfile
    else:
        return None

def get_splicejunction_file(out_dir, data):
    """
    locate the splicejunction file starting from the alignment directory
    """
    samplename = dd.get_sample_name(data)
    sjfile = os.path.join(out_dir, os.pardir, f"{samplename}SJ.out.tab")
    if file_exists(sjfile):
        return sjfile
    else:
        return None

def junction2bed(junction_file):
    """
    reformat the STAR junction file to BED3 format, one end of the splice junction per line
    """
    base, _ = os.path.splitext(junction_file)
    out_file = base + "-minimized.bed"
    if file_exists(out_file):
        return out_file
    if not file_exists(junction_file):
        return None
    with file_transaction(out_file) as tx_out_file:
        with open(junction_file) as in_handle:
            with open(tx_out_file, "w") as out_handle:
                 for line in in_handle:
                    tokens = line.split()
                    chrom, sj1, sj2 = tokens[0:3]
                    if int(sj1) > int(sj2):
                        tmp = sj1
                        sj1 = sj2
                        sj2 = tmp
                    out_handle.write("\t".join([chrom, sj1, sj1]) + "\n")
                    out_handle.write("\t".join([chrom, sj2, sj2]) + "\n")
        minimize = bed.minimize(tx_out_file)
        minimize.saveas(tx_out_file)
    return out_file
