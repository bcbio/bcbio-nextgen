"""Structural variation detection for split and paired reads using lumpy.

Uses speedseq for lumpy integration and samblaster for read preparation:
https://github.com/cc2qe/speedseq
https://github.com/GregoryFaust/samblaster
https://github.com/arq5x/lumpy-sv
"""
import operator
import os
import sys

import toolz as tz

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.ngsalign import postalign
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.structural import delly
from bcbio.variation import vcfutils

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

# Map from numbers used by speedseq to indicate paired and split read evidence
SUPPORT_NUMS = {"1": "PE", "0": "SR"}

def _run_lumpy(full_bams, sr_bams, disc_bams, work_dir, items):
    """Run lumpy-sv, using speedseq pipeline.
    """
    out_file = os.path.join(work_dir, "%s-svs.bedpe"
                            % os.path.splitext(os.path.basename(items[0]["align_bam"]))[0])
    sv_exclude_bed = delly.prepare_exclude_file(items, out_file)
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            with utils.curdir_tmpdir(items[0]) as tmpdir:
                out_base = utils.splitext_plus(tx_out_file)[0]
                full_bams = ",".join(full_bams)
                sr_bams = ",".join(sr_bams)
                disc_bams = ",".join(disc_bams)
                exclude = "-x %s" % sv_exclude_bed if sv_exclude_bed else ""
                cmd = ("speedseq lumpy -v -B {full_bams} -S {sr_bams} -D {disc_bams} {exclude} "
                       "-T {tmpdir} -o {out_base}")
                do.run(cmd.format(**locals()), "speedseq lumpy", items[0])
    return out_file, sv_exclude_bed

def _get_support(parts):
    """Retrieve supporting information for potentially multiple samples.

    Convert speedseqs numbering scheme back into sample and support information.
    sample_ids are generated like 20 or 21, where the first number is sample number
    and the second is the type of supporting evidence.
    """
    out = {}
    for sample_id, read_count in (x.split(",") for x in parts[11].split(":")[-1].split(";")):
        support_type = SUPPORT_NUMS[sample_id[-1]]
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
                        #evidence = set(reduce(operator.add, [x.keys() for x in support.values()]))
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

def _write_samples_to_ids(base_file, items):
    """Write BED file mapping samples to IDs used in the lumpy bedpe output.
    """
    out_file = "%s-samples.bed" % utils.splitext_plus(base_file)[0]
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                for i, data in enumerate(items):
                    sample = tz.get_in(["rgnames", "sample"], data)
                    for sid, stype in SUPPORT_NUMS.items():
                        sample_id = "%s%s" % (i + 1, sid)
                        out_handle.write("%s\t%s\t%s\n" % (sample, sample_id, stype))
    return out_file

def _bedpe_to_vcf(bedpe_file, sconfig_file, items):
    """Convert BEDPE output into a VCF file.
    """
    tovcf_script = do.find_cmd("bedpeToVcf")
    if tovcf_script:
        out_file = "%s.vcf.gz" % utils.splitext_plus(bedpe_file)[0]
        out_nogzip = out_file.replace(".vcf.gz", ".vcf")
        raw_file = "%s-raw.vcf" % utils.splitext_plus(bedpe_file)[0]
        if not utils.file_exists(out_file):
            if not utils.file_exists(raw_file):
                with file_transaction(raw_file) as tx_raw_file:
                    ref_file = tz.get_in(["reference", "fasta", "base"], items[0])
                    cmd = [sys.executable, tovcf_script, "-c", sconfig_file, "-f", ref_file,
                           "-b", bedpe_file, "-o", tx_raw_file]
                    do.run(cmd, "Convert lumpy bedpe output to VCF")
            prep_file = vcfutils.sort_by_ref(raw_file, items[0])
            if not utils.file_exists(out_nogzip):
                utils.symlink_plus(prep_file, out_nogzip)
        out_file = vcfutils.bgzip_and_index(out_nogzip, items[0]["config"])
        return out_file

def _filter_by_bedpe(vcf_file, bedpe_file, data):
    """Add filters to VCF based on pre-filtered bedpe file.
    """
    out_file = "%s-filter%s" % utils.splitext_plus(vcf_file)
    nogzip_out_file = out_file.replace(".vcf.gz", ".vcf")
    if not utils.file_exists(out_file):
        filters = {}
        with open(bedpe_file) as in_handle:
            for line in in_handle:
                parts = line.split("\t")
                name = parts[6]
                cur_filter = parts[-1].strip()
                if cur_filter != "PASS":
                    filters[name] = cur_filter
        with file_transaction(nogzip_out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                with utils.open_gzipsafe(vcf_file) as in_handle:
                    for line in in_handle:
                        if not line.startswith("#"):
                            parts = line.split("\t")
                            cur_id = parts[2].split("_")[0]
                            cur_filter = filters.get(cur_id, "PASS")
                            if cur_filter != "PASS":
                                parts[6] = cur_filter
                            line = "\t".join(parts)
                        out_handle.write(line)
        if out_file.endswith(".gz"):
            vcfutils.bgzip_and_index(nogzip_out_file, data["config"])
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
    pebed_file, exclude_file = _run_lumpy(full_bams, sr_bams, disc_bams, work_dir, items)
    out = []
    sample_config_file = _write_samples_to_ids(pebed_file, items)
    lumpy_vcf = _bedpe_to_vcf(pebed_file, sample_config_file, items)
    for i, data in enumerate(items):
        if "sv" not in data:
            data["sv"] = []
        sample = tz.get_in(["rgnames", "sample"], data)
        sample_bedpe = _filter_by_support(_subset_to_sample(pebed_file, i, data), i)
        if lumpy_vcf:
            sample_vcf = utils.append_stem(lumpy_vcf, "-%s" % sample)
            sample_vcf = _filter_by_bedpe(vcfutils.select_sample(lumpy_vcf, sample, sample_vcf, data["config"]),
                                          sample_bedpe, data)
        else:
            sample_vcf = None
        data["sv"].append({"variantcaller": "lumpy",
                           "vrn_file": sample_vcf,
                           "exclude_file": exclude_file,
                           "bedpe_file": sample_bedpe,
                           "sample_bed": sample_config_file})
        out.append(data)
    return out
