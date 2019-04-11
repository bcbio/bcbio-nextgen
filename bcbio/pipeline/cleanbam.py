"""Clean an input BAM file to work with downstream pipelines.

GATK and Picard based pipelines have specific requirements for
chromosome order, run group information and other BAM formatting.
This provides a pipeline to prepare and resort an input.
"""
import os
import sys

import pysam

from bcbio import bam, broad, utils
from bcbio.bam import ref
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.heterogeneity import chromhacks
from bcbio.ngsalign import novoalign
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do

def fixrg(in_bam, names, ref_file, dirs, data):
    """Fix read group in a file, using samtools addreplacerg.

    addreplacerg does not remove the old read group, causing confusion when
    checking. We use reheader to work around this
    """
    work_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "bamclean", dd.get_sample_name(data)))
    out_file = os.path.join(work_dir, "%s-fixrg.bam" % utils.splitext_plus(os.path.basename(in_bam))[0])
    if not utils.file_exists(out_file):
        out_file = os.path.join(work_dir, "%s-fixrg.bam" % dd.get_sample_name(data))
    if not utils.file_uptodate(out_file, in_bam):
        with file_transaction(data, out_file) as tx_out_file:
            rg_info = novoalign.get_rg_info(names)
            new_header = "%s-header.txt" % os.path.splitext(out_file)[0]
            cores = dd.get_cores(data)
            do.run("samtools view -H {in_bam} | grep -v ^@RG > {new_header}".format(**locals()),
                   "Create empty RG header: %s" % dd.get_sample_name(data))
            cmd = ("samtools reheader {new_header} {in_bam} | "
                   "samtools addreplacerg -@ {cores} -r '{rg_info}' -m overwrite_all -O bam -o {tx_out_file} -")
            do.run(cmd.format(**locals()), "Fix read groups: %s" % dd.get_sample_name(data))
    return out_file

def remove_extracontigs(in_bam, data):
    """Remove extra contigs (non chr1-22,X,Y) from an input BAM.

    These extra contigs can often be arranged in different ways, causing
    incompatibility issues with GATK and other tools. This also fixes the
    read group header as in fixrg.

    This does not yet handle mapping over 1 -> chr1 issues since this requires
    a ton of search/replace which slows down conversion.
    """
    work_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "bamclean", dd.get_sample_name(data)))
    out_file = os.path.join(work_dir, "%s-noextras.bam" % utils.splitext_plus(os.path.basename(in_bam))[0])
    if not utils.file_exists(out_file):
        out_file = os.path.join(work_dir, "%s-noextras.bam" % dd.get_sample_name(data))
    if not utils.file_uptodate(out_file, in_bam):
        with file_transaction(data, out_file) as tx_out_file:
            target_chroms = _target_chroms_and_header(in_bam, data)
            str_chroms = " ".join(target_chroms)
            rg_info = novoalign.get_rg_info(data["rgnames"])
            bcbio_py = sys.executable
            ref_file = dd.get_ref_file(data)
            local_bam = os.path.join(os.path.dirname(tx_out_file), os.path.basename(in_bam))
            cores = dd.get_cores(data)
            utils.symlink_plus(in_bam, local_bam)
            bam.index(local_bam, data["config"])
            cmd = ("samtools view -@ {cores} -h {local_bam} {str_chroms} | "
                   """{bcbio_py} -c 'from bcbio.pipeline import cleanbam; """
                   """cleanbam.fix_header("{ref_file}")' | """
                   "samtools view -@ {cores} -u - | "
                   "samtools addreplacerg -@ {cores} -r '{rg_info}' -m overwrite_all -O bam -o {tx_out_file} - ")
            do.run(cmd.format(**locals()), "bamprep, remove extra contigs: %s" % dd.get_sample_name(data))
    return out_file

def _target_chroms_and_header(bam_file, data):
    """Get a list of chromosomes to target and new updated ref_file header.

    Could potentially handle remapping from chr1 -> 1 but currently disabled due
    to speed issues.
    """
    special_remaps = {"chrM": "MT", "MT": "chrM"}
    target_chroms = dict([(x.name, i) for i, x in enumerate(ref.file_contigs(dd.get_ref_file(data)))
                          if chromhacks.is_autosomal_or_sex(x.name)])
    out_chroms = []
    with pysam.Samfile(bam_file, "rb") as bamfile:
        for bami, bam_contig in enumerate([c["SN"] for c in bamfile.header["SQ"]]):
            if bam_contig in target_chroms:
                target_chrom = bam_contig
            elif bam_contig in special_remaps and special_remaps[bam_contig] in target_chroms:
                target_chrom = special_remaps[bam_contig]
            elif bam_contig.startswith("chr") and bam_contig.replace("chr", "") in target_chroms:
                target_chrom = bam_contig.replace("chr", "")
            elif "chr%s" % bam_contig in target_chroms:
                target_chrom = "chr%s" % bam_contig
            else:
                target_chrom = None
            # target_chrom == bam_contig ensures we don't try chr1 -> 1 style remapping
            if target_chrom and target_chrom == bam_contig:
                # Order not required if dealing with SAM file header fixing
                #assert bami == target_chroms[target_chrom], \
                #    ("remove_extracontigs: Non-matching order of standard contig: %s %s (%s vs %s)" %
                #     (bam_file, target_chrom, bami, target_chroms[target_chrom]))
                out_chroms.append(target_chrom)
    assert out_chroms, ("remove_extracontigs: Did not find any chromosomes in reference file: %s %s" %
                        (bam_file, target_chroms))
    return out_chroms

def fix_header(ref_file):
    added_ref = False
    for line in sys.stdin:
        # skip current read groups, since adding new
        # skip current contigs since adding new sequence dictionary
        if line.startswith(("@RG", "@SQ")):
            pass
        elif not added_ref and not line.startswith("@"):
            for x in ref.file_contigs(ref_file):
                sys.stdout.write("@SQ\tSN:%s\tLN:%s\n" % (x.name, x.size))
            added_ref = True
        else:
            sys.stdout.write(line)

def picard_prep(in_bam, names, ref_file, dirs, data):
    """Prepare input BAM using Picard and GATK cleaning tools.

    - ReorderSam to reorder file to reference
    - AddOrReplaceReadGroups to add read group information and coordinate sort
    - PrintReads to filters to remove problem records:
    - filterMBQ to remove reads with mismatching bases and base qualities
    """
    runner = broad.runner_from_path("picard", data["config"])
    work_dir = utils.safe_makedir(os.path.join(dirs["work"], "bamclean", names["sample"]))
    runner.run_fn("picard_index_ref", ref_file)
    reorder_bam = os.path.join(work_dir, "%s-reorder.bam" %
                               os.path.splitext(os.path.basename(in_bam))[0])
    if not utils.file_exists(reorder_bam):
        reorder_bam = os.path.join(work_dir, "%s-reorder.bam" % dd.get_sample_name(data))
    reorder_bam = runner.run_fn("picard_reorder", in_bam, ref_file, reorder_bam)
    rg_bam = runner.run_fn("picard_fix_rgs", reorder_bam, names)
    return _filter_bad_reads(rg_bam, ref_file, data)

def _filter_bad_reads(in_bam, ref_file, data):
    """Use GATK filter to remove problem reads which choke GATK and Picard.
    """
    bam.index(in_bam, data["config"])
    out_file = "%s-gatkfilter.bam" % os.path.splitext(in_bam)[0]
    if not utils.file_exists(out_file):
        with tx_tmpdir(data) as tmp_dir:
            with file_transaction(data, out_file) as tx_out_file:
                params = [("FixMisencodedBaseQualityReads"
                           if dd.get_quality_format(data, "").lower() == "illumina"
                           else "PrintReads"),
                          "-R", ref_file,
                          "-I", in_bam,
                          "-O", tx_out_file,
                          "-RF", "MatchingBasesAndQualsReadFilter",
                          "-RF", "SeqIsStoredReadFilter",
                          "-RF", "CigarContainsNoNOperator"]
                jvm_opts = broad.get_gatk_opts(data["config"], tmp_dir)
                do.run(broad.gatk_cmd("gatk", jvm_opts, params), "Filter problem reads")
    bam.index(out_file, data["config"])
    return out_file
