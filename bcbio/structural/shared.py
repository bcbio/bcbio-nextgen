"""Shared functionality useful across multiple structural variant callers.

Handles exclusion regions and preparing discordant regions.
"""
from contextlib import closing
import os

import pysam
import toolz as tz

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.bam import callable
from bcbio.ngsalign import postalign
from bcbio.pipeline import shared, config_utils
from bcbio.provenance import do

# ## Case/control

def find_case_control(items):
    """Find case/control items in a population of multiple samples.
    """
    case_phenotypes = set(["affected"])
    control_phenotypes = set(["unaffected"])
    cases = []
    controls = []
    for data in items:
        phenotype = tz.get_in(["metadata", "phenotype"], data)
        if phenotype in case_phenotypes:
            cases.append(data)
        elif phenotype in control_phenotypes:
            controls.append(data)
        else:
            cases.append(data)
    return cases, controls

# ## Prepare exclusion regions (repeats, telomeres, centromeres)

def _get_sv_exclude_file(items):
    """Retrieve SV file of regions to exclude.
    """
    sv_bed = utils.get_in(items[0], ("genome_resources", "variation", "sv_repeat"))
    if sv_bed and os.path.exists(sv_bed):
        return sv_bed

def _get_variant_regions(items):
    """Retrieve variant regions defined in any of the input items.
    """
    return filter(lambda x: x is not None,
                  [tz.get_in(("config", "algorithm", "variant_regions"), data)
                   for data in items
                   if tz.get_in(["config", "algorithm", "coverage_interval"], data) != "genome"])

def has_variant_regions(items, base_file, chrom=None):
    """Determine if we should process this chromosome: needs variant regions defined.
    """
    if chrom:
        all_vrs = _get_variant_regions(items)
        if len(all_vrs) > 0:
            test = shared.subset_variant_regions(tz.first(all_vrs), chrom, base_file, items)
            if test == chrom:
                return False
    return True

def prepare_exclude_file(items, base_file, chrom=None):
    """Prepare a BED file for exclusion, incorporating variant regions and chromosome.

    Excludes locally repetitive regions (if `remove_lcr` is set) and
    centromere regions, both of which contribute to long run times and
    false positive structural variant calls.
    """
    import pybedtools
    out_file = "%s-exclude.bed" % utils.splitext_plus(base_file)[0]
    all_vrs = _get_variant_regions(items)
    ready_region = (shared.subset_variant_regions(tz.first(all_vrs), chrom, base_file, items)
                    if len(all_vrs) > 0 else chrom)
    with shared.bedtools_tmpdir(items[0]):
        # Get a bedtool for the full region if no variant regions
        if ready_region == chrom:
            want_bedtool = callable.get_ref_bedtool(tz.get_in(["reference", "fasta", "base"], items[0]),
                                                    items[0]["config"], chrom)
            lcr_bed = shared.get_lcr_bed(items)
            if lcr_bed:
                want_bedtool = want_bedtool.subtract(pybedtools.BedTool(lcr_bed))
        else:
            want_bedtool = pybedtools.BedTool(ready_region).saveas()
        sv_exclude_bed = _get_sv_exclude_file(items)
        if sv_exclude_bed and len(want_bedtool) > 0:
            want_bedtool = want_bedtool.subtract(sv_exclude_bed).saveas()
        if not utils.file_exists(out_file) and not utils.file_exists(out_file + ".gz"):
            with file_transaction(items[0], out_file) as tx_out_file:
                full_bedtool = callable.get_ref_bedtool(tz.get_in(["reference", "fasta", "base"], items[0]),
                                                        items[0]["config"])
                if len(want_bedtool) > 0:
                    full_bedtool.subtract(want_bedtool).saveas(tx_out_file)
                else:
                    full_bedtool.saveas(tx_out_file)
    return out_file

def get_sv_chroms(items, exclude_file):
    """Retrieve chromosomes to process on, avoiding extra skipped chromosomes.
    """
    import pybedtools
    exclude_regions = {}
    for region in pybedtools.BedTool(exclude_file):
        if int(region.start) == 0:
            exclude_regions[region.chrom] = int(region.end)
    out = []
    with closing(pysam.Samfile(items[0]["work_bam"], "rb")) as pysam_work_bam:
        for chrom, length in zip(pysam_work_bam.references, pysam_work_bam.lengths):
            exclude_length = exclude_regions.get(chrom, 0)
            if exclude_length < length:
                out.append(chrom)
    return out

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
        with tx_tmpdir(data) as tmpdir:
            with file_transaction(data, sr_file) as tx_sr_file:
                with file_transaction(data, disc_file) as tx_disc_file:
                    with file_transaction(data, dedup_file) as tx_dedup_file:
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

def get_split_discordants(data, work_dir):
    """Retrieve full, split and discordant reads, potentially calculating with samblaster as needed.
    """
    dedup_bam, sr_bam, disc_bam = _find_existing_inputs(data["align_bam"])
    if not dedup_bam:
        dedup_bam, sr_bam, disc_bam = _extract_split_and_discordants(data["align_bam"], work_dir, data)
    return dedup_bam, sr_bam, disc_bam

def get_cur_batch(items):
    """Retrieve name of the batch shared between all items in a group.
    """
    batches = []
    for data in items:
        batch = tz.get_in(["metadata", "batch"], data, [])
        batches.append(set(batch) if isinstance(batch, (list, tuple)) else set([batch]))
    combo_batches = reduce(lambda b1, b2: b1.intersection(b2), batches)
    if len(combo_batches) == 1:
        return combo_batches.pop()
    elif len(combo_batches) == 0:
        return None
    else:
        raise ValueError("Found multiple overlapping batches: %s -- %s" % (combo_batches, batches))

def outname_from_inputs(in_files):
    base = os.path.commonprefix(in_files)
    if base.endswith("chr"):
        base = base[:-3]
    while base.endswith(("-", "_", ".")):
        base = base[:-1]
    return base
