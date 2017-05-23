"""Clean an input BAM file to work with downstream pipelines.

GATK and Picard based pipelines have specific requirements for
chromosome order, run group information and other BAM formatting.
This provides a pipeline to prepare and resort an input.
"""
import os
import sys

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
    if not utils.file_uptodate(out_file, in_bam):
        with file_transaction(data, out_file) as tx_out_file:
            rg_info = novoalign.get_rg_info(names)
            new_header = "%s-header.txt" % os.path.splitext(out_file)[0]
            do.run("samtools view -H {in_bam} | grep -v ^@RG > {new_header}".format(**locals()),
                   "Create empty RG header: %s" % dd.get_sample_name(data))
            cmd = ("samtools reheader {new_header} {in_bam} | "
                   "samtools addreplacerg -r '{rg_info}' -m overwrite_all -O bam -o {tx_out_file} -")
            do.run(cmd.format(**locals()), "Fix read groups: %s" % dd.get_sample_name(data))
    return out_file

def remove_extracontigs(in_bam, data):
    """Remove extra contigs (non chr1-22,X,Y) from an input BAM.

    These extra contigs can often be arranged in different ways, causing
    incompatibility issues with GATK and other tools. This also fixes the
    read group header as in fixrg.
    """
    work_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "bamclean", dd.get_sample_name(data)))
    out_file = os.path.join(work_dir, "%s-noextras.bam" % utils.splitext_plus(os.path.basename(in_bam))[0])
    if not utils.file_uptodate(out_file, in_bam):
        with file_transaction(data, out_file) as tx_out_file:
            target_chroms = [x.name for x in ref.file_contigs(dd.get_ref_file(data))
                             if chromhacks.is_autosomal_or_sex(x.name)]
            str_chroms = " ".join(target_chroms)
            comma_chroms = ",".join(target_chroms)
            rg_info = novoalign.get_rg_info(data["rgnames"])
            bcbio_py = sys.executable
            cmd = ("samtools view -h {in_bam} {str_chroms} | "
                   """{bcbio_py} -c 'from bcbio.pipeline import cleanbam; """
                   """cleanbam.fix_header("{comma_chroms}")' | """
                   "samtools view -u - | "
                   "samtools addreplacerg -r '{rg_info}' -m overwrite_all -O bam -o {tx_out_file} - ")
            do.run(cmd.format(**locals()), "bamprep, remove extra contigs: %s" % dd.get_sample_name(data))
    return out_file

def fix_header(chroms):
    target_chroms = chroms.split(",")
    for line in sys.stdin:
        if line.startswith("@RG"):
            pass  # skip current read groups, since adding new
        elif line.startswith("@SQ"):
            contig = [x for x in line.split("\t") if x.startswith("SN:")][0].split(":")[-1]
            if contig in target_chroms:
                sys.stdout.write(line)
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
                params = ["-T", "PrintReads",
                          "-R", ref_file,
                          "-I", in_bam,
                          "--out", tx_out_file,
                          "--filter_mismatching_base_and_quals",
                          "--filter_bases_not_stored",
                          "--filter_reads_with_N_cigar"]
                if dd.get_quality_format(data, "").lower() == "illumina":
                    params.append("--fix_misencoded_quality_scores")
                jvm_opts = broad.get_gatk_framework_opts(data["config"], tmp_dir)
                do.run(broad.gatk_cmd("gatk-framework", jvm_opts, params), "Filter problem reads")
    bam.index(out_file, data["config"])
    return out_file
