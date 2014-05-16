"""Structural variation detection for split and paired reads using lumpy.

Uses speedseq for lumpy integration and samblaster for read preparation:
https://github.com/cc2qe/speedseq
https://github.com/GregoryFaust/samblaster
https://github.com/arq5x/lumpy-sv
"""
import operator
import os

import toolz as tz

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.ngsalign import postalign
from bcbio.pipeline import config_utils, shared
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
                                     3, "decrease").upper()
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
                sv_exclude_bed = delly.prepare_exclude_file(items, out_file)
                exclude = "-x %s" % sv_exclude_bed if sv_exclude_bed else ""
                cmd = ("speedseq lumpy -v -B {full_bams} -S {sr_bams} -D {disc_bams} {exclude} "
                       "-T {tmpdir} -o {out_base}")
                do.run(cmd.format(**locals()), "speedseq lumpy", items[0])
    return out_file

def _get_support(parts):
    """Retrieve supporting information for potentially multiple samples.

    Convert speedseqs numbering scheme back into sample and support information.
    sample_ids are generated like 20 or 21, where the first number is sample number
    and the second is the type of supporting evidence.
    """
    support_nums = {"1": "discordant", "0": "split-reads"}
    out = {}
    for sample_id, read_count in (x.split(",") for x in parts[11].split(":")[-1].split(";")):
        support_type = support_nums[sample_id[-1]]
        sample_id = int(sample_id[:-1]) - 1
        out = tz.update_in(out, [sample_id, support_type], lambda x: x + int(read_count), 0)
    return out

def _subset_to_sample(orig_file, index, data):
    """Subset population based calls to those supported within a single sample.
    """
    out_file = utils.append_stem(orig_file, "-" + data["rgnames"]["sample"])
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            with open(orig_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    for parts in (l.rstrip().split("\t") for l in in_handle):
                        support = _get_support(parts)
                        if index in support:
                            out_handle.write("\t".join(parts) + "\n")
    return out_file

def _filter_by_support(orig_file, index):
    """Filter call file based on supporting evidence, adding pass/filter annotations to BEDPE.

    Filters based on the following criteria:
      - Minimum read support for the call.
    Other filters not currently applied due to being too restrictive:
      - Multiple forms of evidence in any sample (split and paired end)
    """
    min_read_count = 4
    out_file = utils.append_stem(orig_file, "-filter")
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            with open(orig_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    for parts in (l.rstrip().split("\t") for l in in_handle):
                        support = _get_support(parts)
                        evidence = set(reduce(operator.add, [x.keys() for x in support.values()]))
                        read_count = reduce(operator.add, support[index].values())
                        if read_count < min_read_count:
                            lfilter = "ReadCountSupport"
                        #elif len(evidence) < 2:
                        #    lfilter = "ApproachSupport"
                        else:
                            lfilter = "PASS"
                        parts.append(lfilter)
                        out_handle.write("\t".join(parts) + "\n")
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
    for i, data in enumerate(items):
        if "sv" not in data:
            data["sv"] = []
        data["sv"].append({"variantcaller": "lumpy",
                           "vrn_file": _filter_by_support(_subset_to_sample(pebed_file, i, data), i)})
        out.append(data)
    return out
