"""Structural variation detection for split and paired reads using lumpy.

Uses smoove for automating lumpy variant calling:
https://github.com/brentp/smoove
https://github.com/arq5x/lumpy-sv
"""
import collections
import os
import re
import subprocess

import six
import vcf

from bcbio import utils
from bcbio.bam import ref
from bcbio.log import logger
from bcbio.distributed.transaction import file_transaction
from bcbio.heterogeneity import chromhacks
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.structural import shared as sshared
from bcbio.variation import effects, vcfutils, vfilter

# ## Lumpy main

def _run_smoove(full_bams, sr_bams, disc_bams, work_dir, items):
    """Run lumpy-sv using smoove.
    """
    batch = sshared.get_cur_batch(items)
    ext = "-%s-svs" % batch if batch else "-svs"
    name = "%s%s" % (dd.get_sample_name(items[0]), ext)
    out_file = os.path.join(work_dir, "%s-smoove.genotyped.vcf.gz" % name)
    sv_exclude_bed = sshared.prepare_exclude_file(items, out_file)
    old_out_file = os.path.join(work_dir, "%s%s-prep.vcf.gz"
                                % (os.path.splitext(os.path.basename(items[0]["align_bam"]))[0], ext))
    if utils.file_exists(old_out_file):
        return old_out_file, sv_exclude_bed
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            cores = dd.get_num_cores(items[0])
            out_dir = os.path.dirname(tx_out_file)
            ref_file = dd.get_ref_file(items[0])
            full_bams = " ".join(_prepare_smoove_bams(full_bams, sr_bams, disc_bams, items,
                                                      os.path.dirname(tx_out_file)))
            std_excludes = ["~^GL", "~^HLA", "~_random", "~^chrUn", "~alt", "~decoy"]
            def _is_std_exclude(n):
                clean_excludes = [x.replace("~", "").replace("^", "") for x in std_excludes]
                return any([n.startswith(x) or n.endswith(x) for x in clean_excludes])
            exclude_chrs = [c.name for c in ref.file_contigs(ref_file)
                            if not chromhacks.is_nonalt(c.name) and not _is_std_exclude(c.name)]
            exclude_chrs = "--excludechroms '%s'" % ",".join(std_excludes + exclude_chrs)
            exclude_bed = ("--exclude %s" % sv_exclude_bed) if utils.file_exists(sv_exclude_bed) else ""
            tempdir = os.path.dirname(tx_out_file)
            cmd = ("export TMPDIR={tempdir} && "
                   "smoove call --processes {cores} --genotype --removepr --fasta {ref_file} "
                   "--name {name} --outdir {out_dir} "
                   "{exclude_bed} {exclude_chrs} {full_bams}")
            with utils.chdir(tempdir):
                try:
                    do.run(cmd.format(**locals()), "smoove lumpy calling", items[0])
                except subprocess.CalledProcessError as msg:
                    if _allowed_errors(msg):
                        vcfutils.write_empty_vcf(tx_out_file, config=items[0]["config"],
                                                 samples=[dd.get_sample_name(d) for d in items])
                    else:
                        logger.exception()
                        raise
    vcfutils.bgzip_and_index(out_file, items[0]["config"])
    return out_file, sv_exclude_bed

def _prepare_smoove_bams(full_bams, sr_bams, disc_bams, items, tx_work_dir):
    """Prepare BAMs for smoove, linking in pre-existing split/disc BAMs if present.

    Smoove can use pre-existing discordant and split BAMs prepared by samblaster
    if present as $sample.split.bam and $sample.disc.bam.
    """
    input_dir = utils.safe_makedir(tx_work_dir)
    out = []
    for full, sr, disc, data in zip(full_bams, sr_bams, disc_bams, items):
        if sr and disc:
            new_full = os.path.join(input_dir, "%s.bam" % dd.get_sample_name(data))
            new_sr = os.path.join(input_dir, "%s.split.bam" % dd.get_sample_name(data))
            new_disc = os.path.join(input_dir, "%s.disc.bam" % dd.get_sample_name(data))
            utils.symlink_plus(full, new_full)
            utils.symlink_plus(sr, new_sr)
            utils.symlink_plus(disc, new_disc)
            out.append(new_full)
        else:
            out.append(full)
    return out

def _allowed_errors(msg):
    if six.PY3:
        msg = str(msg)
    else:
        msg = unicode(msg).encode("ascii", "replace")
    allowed = ["covmed: not enough reads to sample for bam stats",
               "missing pair end parameters:",
               "mean stdev read_length min_non_overlap"]
    return any([len(re.findall(m, msg)) > 0 for m in allowed])

def _filter_by_support(in_file, data):
    """Filter call file based on supporting evidence, adding FILTER annotations to VCF.

    Filters based on the following criteria:
      - Minimum read support for the call (SU = total support)
      - Large calls need split read evidence.
    """
    rc_filter = ("FORMAT/SU < 4 || "
                 "(FORMAT/SR == 0 && FORMAT/SU < 15 && ABS(SVLEN)>50000) || "
                 "(FORMAT/SR == 0 && FORMAT/SU < 5 && ABS(SVLEN)<2000) || "
                 "(FORMAT/SR == 0 && FORMAT/SU < 15 && ABS(SVLEN)<300)")
    return vfilter.cutoff_w_expression(in_file, rc_filter, data, name="ReadCountSupport",
                                       limit_regions=None)

def _filter_by_background(base_name, back_samples, gt_vcfs, data):
    """Filter base samples, marking any also present in the background.
    """
    filtname = "InBackground"
    filtdoc = "Variant also present in background samples with same genotype"
    orig_vcf = gt_vcfs[base_name]
    out_file = "%s-backfilter.vcf" % (utils.splitext_plus(orig_vcf)[0])
    if not utils.file_exists(out_file) and not utils.file_exists(out_file + ".gz"):
        with file_transaction(data, out_file) as tx_out_file:
            with utils.open_gzipsafe(orig_vcf) as in_handle:
                inp = vcf.Reader(in_handle, orig_vcf)
                inp.filters[filtname] = vcf.parser._Filter(filtname, filtdoc)
                with open(tx_out_file, "w") as out_handle:
                    outp = vcf.Writer(out_handle, inp)
                    for rec in inp:
                        if _genotype_in_background(rec, base_name, back_samples):
                            rec.add_filter(filtname)
                        outp.write_record(rec)
    if utils.file_exists(out_file + ".gz"):
        out_file = out_file + ".gz"
    gt_vcfs[base_name] = vcfutils.bgzip_and_index(out_file, data["config"])
    return gt_vcfs

def _genotype_in_background(rec, base_name, back_samples):
    """Check if the genotype in the record of interest is present in the background records.
    """
    def passes(rec):
        return not rec.FILTER or len(rec.FILTER) == 0
    return (passes(rec) and
            any(rec.genotype(base_name).gt_alleles == rec.genotype(back_name).gt_alleles
                for back_name in back_samples))

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "lumpy"))

def run(items):
    """Perform detection of structural variations with lumpy.
    """
    paired = vcfutils.get_paired(items)
    work_dir = _sv_workdir(paired.tumor_data if paired and paired.tumor_data else items[0])
    previous_evidence = {}
    full_bams, sr_bams, disc_bams = [], [], []
    for data in items:
        full_bams.append(dd.get_align_bam(data))
        sr_bam, disc_bam = sshared.find_existing_split_discordants(data)
        sr_bams.append(sr_bam)
        disc_bams.append(disc_bam)
        cur_dels, cur_dups = _bedpes_from_cnv_caller(data, work_dir)
        previous_evidence[dd.get_sample_name(data)] = {}
        if cur_dels and utils.file_exists(cur_dels):
            previous_evidence[dd.get_sample_name(data)]["dels"] = cur_dels
        if cur_dups and utils.file_exists(cur_dups):
            previous_evidence[dd.get_sample_name(data)]["dups"] = cur_dups
    lumpy_vcf, exclude_file = _run_smoove(full_bams, sr_bams, disc_bams, work_dir, items)
    lumpy_vcf = sshared.annotate_with_depth(lumpy_vcf, items)
    gt_vcfs = {}
    # Retain paired samples with tumor/normal genotyped in one file
    if paired and paired.normal_name:
        batches = [[paired.tumor_data, paired.normal_data]]
    else:
        batches = [[x] for x in items]

    for batch_items in batches:
        for data in batch_items:
            gt_vcfs[dd.get_sample_name(data)] = _filter_by_support(lumpy_vcf, data)
    if paired and paired.normal_name:
        gt_vcfs = _filter_by_background(paired.tumor_name, [paired.normal_name], gt_vcfs, paired.tumor_data)
    out = []
    upload_counts = collections.defaultdict(int)
    for data in items:
        if "sv" not in data:
            data["sv"] = []
        vcf_file = gt_vcfs.get(dd.get_sample_name(data))
        if vcf_file:
            effects_vcf, _ = effects.add_to_vcf(vcf_file, data, "snpeff")
            data["sv"].append({"variantcaller": "lumpy",
                               "vrn_file": effects_vcf or vcf_file,
                               "do_upload": upload_counts[vcf_file] == 0,  # only upload a single file per batch
                               "exclude_file": exclude_file})
            upload_counts[vcf_file] += 1
        out.append(data)
    return out

def _bedpes_from_cnv_caller(data, work_dir):
    """Retrieve BEDPEs deletion and duplications from CNV callers.

    Currently integrates with CNVkit.
    """
    supported = set(["cnvkit"])
    cns_file = None
    for sv in data.get("sv", []):
        if sv["variantcaller"] in supported and "cns" in sv and "lumpy_usecnv" in dd.get_tools_on(data):
            cns_file = sv["cns"]
            break
    if not cns_file:
        return None, None
    else:
        out_base = os.path.join(work_dir, utils.splitext_plus(os.path.basename(cns_file))[0])
        out_dels = out_base + "-dels.bedpe"
        out_dups = out_base + "-dups.bedpe"
        if not os.path.exists(out_dels) or not os.path.exists(out_dups):
            with file_transaction(data, out_dels, out_dups) as (tx_out_dels, tx_out_dups):
                try:
                    cnvanator_path = config_utils.get_program("cnvanator_to_bedpes.py", data)
                except config_utils.CmdNotFound:
                    return None, None
                cmd = [cnvanator_path, "-c", cns_file, "--cnvkit",
                        "--del_o=%s" % tx_out_dels, "--dup_o=%s" % tx_out_dups,
                        "-b", "250"]  # XXX Uses default piece size for CNVkit. Right approach?
                do.run(cmd, "Prepare CNVkit as input for lumpy", data)
        return out_dels, out_dups
