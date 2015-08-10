import os
import sys
import shutil
import subprocess
import contextlib
from distutils.version import LooseVersion

from bcbio.pipeline import config_utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.utils import (safe_makedir, file_exists, is_gzipped)
from bcbio.provenance import do
from bcbio import bam, utils
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.ngsalign import postalign

CLEANUP_FILES = ["Aligned.out.sam", "Log.out", "Log.progress.out"]
ALIGN_TAGS = ["NH", "HI", "NM", "MD", "AS"]

def align(fastq_file, pair_file, ref_file, names, align_dir, data):
    max_hits = 10
    srna = True if data["analysis"].lower().startswith("smallrna-seq") else False
    srna_opts = ""
    if srna:
        max_hits = 1000
        srna_opts = "--alignIntronMax 1"
    config = data["config"]
    out_prefix = os.path.join(align_dir, dd.get_lane(data))
    out_file = out_prefix + "Aligned.out.sam"
    out_dir = os.path.join(align_dir, "%s_star" % dd.get_lane(data))

    if not ref_file:
        logger.error("STAR index not found. We don't provide the STAR indexes "
                     "by default because they are very large. You can install "
                     "the index for your genome with: bcbio_nextgen.py upgrade "
                     "--aligners star --genomes genome-build-name --data")
        sys.exit(1)

    final_out = os.path.join(out_dir, "{0}.bam".format(names["sample"]))
    if file_exists(final_out):
        data = _update_data(final_out, out_dir, names, data)
        return data

    star_path = config_utils.get_program("STAR", config)
    fastq = " ".join([fastq_file, pair_file]) if pair_file else fastq_file
    num_cores = config["algorithm"].get("num_cores", 1)

    safe_makedir(align_dir)
    cmd = ("{star_path} --genomeDir {ref_file} --readFilesIn {fastq} "
           "--runThreadN {num_cores} --outFileNamePrefix {out_prefix} "
           "--outReadsUnmapped Fastx --outFilterMultimapNmax {max_hits} "
           "--outStd SAM {srna_opts} "
           "--outSAMunmapped Within --outSAMattributes %s" % " ".join(ALIGN_TAGS))
    cmd = cmd + " --readFilesCommand zcat " if is_gzipped(fastq_file) else cmd
    cmd += _read_group_option(names)
    fusion_mode = utils.get_in(data, ("config", "algorithm", "fusion_mode"), False)
    if fusion_mode:
        cmd += " --chimSegmentMin 15 --chimJunctionOverhangMin 15"
    strandedness = utils.get_in(data, ("config", "algorithm", "strandedness"),
                                "unstranded").lower()
    if strandedness == "unstranded" and not srna:
        cmd += " --outSAMstrandField intronMotif "

    if dd.get_transcriptome_align(data) and not is_transcriptome_broken():
        cmd += " --quantMode TranscriptomeSAM "

    with file_transaction(data, final_out) as tx_final_out:
        cmd += " | " + postalign.sam_to_sortbam_cl(data, tx_final_out)
        run_message = "Running STAR aligner on %s and %s" % (fastq_file, ref_file)
        do.run(cmd.format(**locals()), run_message, None)

    data = _update_data(final_out, out_dir, names, data)
    return data

def _update_data(align_file, out_dir, names, data):
    data = dd.set_work_bam(data, align_file)
    data = dd.set_align_bam(data, align_file)
    if dd.get_transcriptome_align(data) and not is_transcriptome_broken():
        transcriptome_file = _move_transcriptome_file(out_dir, names)
        data = dd.set_transcriptome_bam(data, transcriptome_file)
    return data

def _move_transcriptome_file(out_dir, names):
    out_file = os.path.join(out_dir, "{0}.transcriptome.bam".format(names["sample"]))
    if not file_exists(out_file):
        tmp_file = os.path.join(out_dir, os.pardir,
                                "{0}Aligned.toTranscriptome.out.bam".format(names["lane"]))
        shutil.move(tmp_file, out_file)
    return out_file

def _read_group_option(names):
    rg_id = names["rg"]
    rg_sample = names["sample"]
    rg_library = names["pl"]
    rg_platform_unit = names["pu"]

    return (" --outSAMattrRGline ID:{rg_id} PL:{rg_library} "
            "PU:{rg_platform_unit} SM:{rg_sample} ").format(**locals())

def _get_quality_format(config):
    qual_format = config["algorithm"].get("quality_format", None)
    if qual_format.lower() == "illumina":
        return "fastq-illumina"
    elif qual_format.lower() == "solexa":
        return "fastq-solexa"
    else:
        return "fastq-sanger"

def remap_index_fn(ref_file):
    """Map sequence references to equivalent star indexes
    """
    return os.path.join(os.path.dirname(os.path.dirname(ref_file)), "star")

def index(ref_file, out_dir, data):
    """Create a STAR index in the defined reference directory.
    """
    (ref_dir, local_file) = os.path.split(ref_file)
    gtf_file = os.path.join(ref_dir, os.pardir, "rnaseq", "ref-transcripts.gtf")
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

def get_star_version():
    cmd = "STAR --version"
    subp = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            shell=True)
    with contextlib.closing(subp.stdout) as stdout:
        for line in stdout:
            if "STAR_" in line:
                version = line.split("STAR_")[1].strip()
    return version

def is_transcriptome_broken():
    """
    mapping to the transcriptome causes segfaults, but will be supported later
    until it is fixed, skip this and use the fallback mapping instead
    """
    version = get_star_version()
    return LooseVersion(version) >= LooseVersion("2.4.0g")
