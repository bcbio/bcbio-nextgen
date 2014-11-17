import os
import sys
import shutil

from bcbio.pipeline import config_utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.utils import (safe_makedir, file_exists, is_gzipped)
from bcbio.provenance import do
from bcbio import bam, utils
from bcbio.log import logger
from bcbio.pipeline import datadict as dd

CLEANUP_FILES = ["Aligned.out.sam", "Log.out", "Log.progress.out"]
ALIGN_TAGS = ["NH", "HI", "NM", "MD", "AS"]

def align(fastq_file, pair_file, ref_file, names, align_dir, data):
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
        return final_out
    star_path = config_utils.get_program("STAR", config)
    fastq = " ".join([fastq_file, pair_file]) if pair_file else fastq_file
    num_cores = config["algorithm"].get("num_cores", 1)

    safe_makedir(align_dir)
    cmd = ("{star_path} --genomeDir {ref_file} --readFilesIn {fastq} "
           "--runThreadN {num_cores} --outFileNamePrefix {out_prefix} "
           "--outReadsUnmapped Fastx --outFilterMultimapNmax 10 "
           "--outStd SAM "
           "--outSAMunmapped Within --outSAMattributes %s" % " ".join(ALIGN_TAGS))
    cmd = cmd + " --readFilesCommand zcat " if is_gzipped(fastq_file) else cmd
    cmd += _read_group_option(names)
    fusion_mode = utils.get_in(data, ("config", "algorithm", "fusion_mode"), False)
    if fusion_mode:
        cmd += " --chimSegmentMin 15 --chimJunctionOverhangMin 15"
    strandedness = utils.get_in(data, ("config", "algorithm", "strandedness"),
                                "unstranded").lower()
    if strandedness == "unstranded":
        cmd += " --outSAMstrandField intronMotif "

    if dd.get_rsem(data):
        cmd += " --quantMode TranscriptomeSAM "

    with tx_tmpdir(data) as tmp_dir:
        sam_to_bam = bam.sam_to_bam_stream_cmd(config)
        sort = bam.sort_cmd(config, tmp_dir)
        cmd += "| {sam_to_bam} | {sort} -o {tx_final_out} "
        run_message = "Running STAR aligner on %s and %s" % (fastq_file, ref_file)
        with file_transaction(data, final_out) as tx_final_out:
            do.run(cmd.format(**locals()), run_message, None)

    if dd.get_rsem(data):
        transcriptome_file = _move_transcriptome_file(out_dir, names)
    return final_out

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
