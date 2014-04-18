"""Structural variation detection for split and paired reads using lumpy.

Uses speedseq for lumpy integration and samblaster for read preparation:
https://github.com/cc2qe/speedseq
https://github.com/GregoryFaust/samblaster
https://github.com/arq5x/lumpy-sv
"""
import os

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.ngsalign import postalign
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.structural import delly

# ## Read preparation

def _extract_split_and_discordants(in_bam, work_dir, data):
    """Retrieve split-read alignments from input BAM file.
    """
    dedup_file = os.path.join(work_dir, "%s-dedup.bam" % os.path.splitext(os.path.basename(in_bam))[0])
    sr_file = os.path.join(work_dir, "%s-sr.bam" % os.path.splitext(os.path.basename(in_bam))[0])
    disc_file = os.path.join(work_dir, "%s-disc.bam" % os.path.splitext(os.path.basename(in_bam))[0])
    samtools = config_utils.get_program("samtools", data["config"])
    cores = utils.get_in(data, ("config", "algorithm", "num_cores"), 1)
    resources = config_utils.get_resources("sambamba", data["config"])
    mem = config_utils.adjust_memory(resources.get("memory", "2G"),
                                     3, "decrease")
    if not utils.file_exists(sr_file) or not utils.file_exists(disc_file) or utils.file_exists(dedup_file):
        with utils.curdir_tmpdir(data) as tmpdir:
            with file_transaction(sr_file) as tx_sr_file:
                with file_transaction(disc_file) as tx_disc_file:
                    with file_transaction(dedup_file) as tx_dedup_file:
                        samblaster_cl = postalign.samblaster_dedup_sort(data, tmpdir, tx_dedup_file,
                                                                        tx_sr_file, tx_disc_file)
                        out_base = os.path.join(tmpdir, "%s-namesort" % os.path.splitext(in_bam)[0])
                        cmd = ("{samtools} sort -n -o -@ {cores} -m {mem} {in_bam} {out_base} | "
                               "{samtools} view -h - | ")
                        cmd = cmd.format(**locals()) + samblaster_cl
                        do.run(cmd, "samblaster: split and discordant reads", data)
    for fname in [sr_file, disc_file, dedup_file]:
        bam.index(fname, data["config"])
    return dedup_file, sr_file, disc_file

def _find_existing_inputs(in_bam):
    """Check for pre-calculated split reads and discordants done as part of alignment streaming.
    """
    sr_file = "%s-sr.bam" % os.path.splitext(in_bam)[0]
    disc_file = "%s-disc.bam" % os.path.splitext(in_bam)[0]
    if utils.file_exists(sr_file) and utils.file_exists(disc_file):
        return in_bam, sr_file, disc_file
    else:
        return None, None, None

# ## Lumpy main

def _run_lumpy(full_bams, sr_bams, disc_bams, work_dir, items):
    """Run lumpy-sv, using speedseq pipeline.
    """
    out_file = os.path.join(work_dir, "%s-svs.bedpe"
                            % os.path.splitext(os.path.basename(items[0]["align_bam"]))[0])
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            with utils.curdir_tmpdir(items[0]) as tmpdir:
                out_base = utils.splitext_plus(tx_out_file)[0]
                full_bams = ",".join(full_bams)
                sr_bams = ",".join(sr_bams)
                disc_bams = ",".join(disc_bams)
                sv_exclude_bed = delly.get_sv_exclude_file(items)
                exclude = "-x %s" % sv_exclude_bed if sv_exclude_bed else ""
                cmd = ("speedseq lumpy -B {full_bams} -S {sr_bams} -D {disc_bams} {exclude} "
                       "-T {tmpdir} -o {out_base}")
                do.run(cmd.format(**locals()), "speedseq lumpy", items[0])
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
        dedup_bam, sr_bam, disc_bam = _find_existing_inputs(data["align_bam"])
        if not dedup_bam:
            dedup_bam, sr_bam, disc_bam = _extract_split_and_discordants(data["align_bam"], work_dir, data)
        full_bams.append(dedup_bam)
        sr_bams.append(sr_bam)
        disc_bams.append(disc_bam)
    pebed_file = _run_lumpy(full_bams, sr_bams, disc_bams, work_dir, items)
    out = []
    for data in items:
        if "sv" not in data:
            data["sv"] = {}
        data["sv"]["lumpy"] = pebed_file
        out.append(data)
    return out
