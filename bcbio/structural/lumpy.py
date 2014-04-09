"""Structural variation detection for split and paired reads using lumpy.

Uses speedseq for lumpy integration and samblaster for read preparation:
https://github.com/cc2qe/speedseq
https://github.com/GregoryFaust/samblaster
https://github.com/arq5x/lumpy-sv
"""
import os

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.provenance import do

# ## Read preparation

def _extract_split_and_discordants(in_bam, work_dir, data):
    """Retrieve split-read alignments from input BAM file.
    """
    sr_file = os.path.join(work_dir, "%s-sr.bam" % os.path.splitext(os.path.basename(in_bam))[0])
    disc_file = os.path.join(work_dir, "%s-disc.bam" % os.path.splitext(os.path.basename(in_bam))[0])
    samblaster = config_utils.get_program("samblaster", data["config"])
    sambamba = config_utils.get_program("sambamba", data["config"])
    cores = utils.get_in(data, ("config", "algorithm", "num_cores"), 1)
    resources = config_utils.get_resources("sambamba", data["config"])
    mem = config_utils.adjust_memory(resources.get("memory", "2G"),
                                     3, "decrease")
    if not utils.file_exists(sr_file) or not utils.file_exists(disc_file):
        with file_transaction(sr_file) as tx_sr_file:
            with file_transaction(disc_file) as tx_disc_file:
                with utils.curdir_tmpdir() as tmpdir:
                    tobam_cmd = ("{sambamba} view -S -f bam -l 0 /dev/stdin | "
                                 "{sambamba} sort -t {cores} -m {mem} --tmpdir {tmpdir} "
                                 "-o {out_file} /dev/stdin")
                    splitter_cmd = tobam_cmd.format(out_file=tx_sr_file, **locals())
                    discordant_cmd = tobam_cmd.format(out_file=tx_disc_file, **locals())
                    cmd = ("{sambamba} sort -t {cores} -m {mem} --tmpdir={tmpdir} "
                           "-n -o /dev/stdout -l 0 {in_bam} | "
                           "{sambamba} view -h /dev/stdin | "
                           "{samblaster} --splitterFile >({splitter_cmd}) --discordantFile >({discordant_cmd}) "
                           "-o /dev/null")
                    do.run(cmd.format(**locals()), "samblaster: split and discordant reads", data)
    return sr_file, disc_file

# ## Lumpy main

def _run_lumpy(full_bams, sr_bams, disc_bams, work_dir, data):
    """Run lumpy-sv, using speedseq pipeline.
    """
    out_file = os.path.join(work_dir, "%s-svs.bedpe"
                            % os.path.splitext(os.path.basename(data["work_bam"]))[0])
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            with utils.curdir_tmpdir() as tmpdir:
                out_base = utils.splitext_plus(tx_out_file)[0]
                full_bams = ",".join(full_bams)
                sr_bams = ",".join(sr_bams)
                disc_bams = ",".join(disc_bams)
                cmd = ("speedseq lumpy -B {full_bams} -S {sr_bams} -D {disc_bams} "
                       "-T {tmpdir} -o {out_base}")
                do.run(cmd.format(**locals()), "speedseq lumpy", data)
    return out_file

def run(items):
    """Perform detection of structural variations with lumpy, using bwa-mem alignment.
    """
    if not all(utils.get_in(data, ("config", "algorithm", "aligner")) == "bwa" for data in items):
        raise ValueError("Require bwa-mem alignment input for lumpy structural variation detection")
    work_dir = utils.safe_makedir(os.path.join(items[0]["dirs"]["work"], "structural", items[0]["name"][-1],
                                               "lumpy"))
    full_bams, sr_bams, disc_bams = [], [], []
    for data in items:
        full_bams.append(data["work_bam"])
        sr_bam, disc_bam = _extract_split_and_discordants(data["work_bam"], work_dir, data)
        sr_bams.append(sr_bam)
        disc_bams.append(disc_bam)
    pebed_file = _run_lumpy(full_bams, sr_bams, disc_bams, work_dir, data)
    out = []
    for data in items:
        if "sv" not in data:
            data["sv"] = {}
        data["sv"]["lumpy"] = pebed_file
        out.append(data)
    return out
