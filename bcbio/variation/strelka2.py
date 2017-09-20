"""Germline and somatic calling with Strelka2: https://github.com/illumina/strelka
"""
import os
import sys

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import shared
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import annotation, bedutils, ploidy, vcfutils

def run(align_bams, items, ref_file, assoc_files, region=None, out_file=None):
    """Run strelka2 variant calling, either paired tumor/normal or germline calling.
    """
    if vcfutils.is_paired_analysis(align_bams, items):
        paired = vcfutils.get_paired_bams(align_bams, items)
        assert paired.normal_bam, "Strelka2 requires a normal sample"
        call_file = _run_somatic(paired, ref_file, assoc_files, region, out_file)
    else:
        call_file = _run_germline(align_bams, items, ref_file,
                                  assoc_files, region, out_file)
    return call_file

def _get_region_bed(region, items, out_file):
    variant_regions = bedutils.merge_overlaps(bedutils.population_variant_regions(items), items[0])
    target = shared.subset_variant_regions(variant_regions, region, out_file, items)
    if not target:
        raise ValueError("Need BED input for strelka2 regions: %s %s" % (region, target))
    if not isinstance(target, basestring) or not os.path.isfile(target):
        chrom, start, end = target
        target = "%s-regions.bed" % utils.splitext_plus(out_file)[0]
        with file_transaction(items[0], target) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                out_handle.write("%s\t%s\t%s\n" % (chrom, start, end))
    return bedutils.merge_overlaps(target, items[0], out_dir=os.path.dirname(out_file)) + ".gz"

def _get_ploidy(region, items, base_file):
    samples = [dd.get_sample_name(d) for d in items]
    ploidies = [ploidy.get_ploidy([d], region) for d in items]
    out_file = "%s-ploidy.vcf" % utils.splitext_plus(base_file)[0]
    if not utils.file_exists(out_file) and not utils.file_exists(out_file + ".gz"):
        with file_transaction(items[0], out_file) as tx_outfile:
            with open(tx_outfile, "w") as h:
                h.write("##fileformat=VCFv4.1\n")
                h.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n')
                h.write('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">\n')
                h.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
                h.write("\t".join([region[0], str(region[1]), ".", "N", "<CNV>", ".", ".",
                                   "END=%s" % region[2], "CN"] + [str(x) for x in ploidies]) + "\n")
    return vcfutils.bgzip_and_index(out_file, items[0]["config"])

def _configure_germline(align_bams, items, ref_file, region, out_file, tx_work_dir):
    utils.safe_makedir(tx_work_dir)
    cmd = [sys.executable, os.path.realpath(utils.which("configureStrelkaGermlineWorkflow.py"))]
    cmd += ["--referenceFasta=%s" % ref_file,
            "--callRegions=%s" % _get_region_bed(region, items, out_file),
            "--ploidy=%s" % _get_ploidy(region, items, out_file),
            "--runDir=%s" % tx_work_dir]
    cmd += ["--bam=%s" % b for b in align_bams]
    if any(dd.get_coverage_interval(d) not in ["genome"] for d in items):
        cmd += ["--targeted"]
    do.run(cmd, "Configure Strelka2 germline calling: %s" % (", ".join([dd.get_sample_name(d) for d in items])))
    return os.path.join(tx_work_dir, "runWorkflow.py")

def _run_germline(align_bams, items, ref_file, assoc_files, region, out_file):
    if not utils.file_exists(out_file):
        work_dir = "%s-work" % utils.splitext_plus(out_file)[0]
        with file_transaction(items[0], work_dir) as tx_work_dir:
            workflow_file = _configure_germline(align_bams, items, ref_file, region, out_file, tx_work_dir)
            _run_workflow(items[0], workflow_file, tx_work_dir)
        raw_file = os.path.join(work_dir, "results", "variants", "variants.vcf.gz")
        out_file = annotation.annotate_nongatk_vcf(raw_file, align_bams, assoc_files.get("dbsnp"),
                                                   ref_file, items[0], out_file)
    return vcfutils.bgzip_and_index(out_file, items[0]["config"])

def _configure_somatic(paired, ref_file, region, out_file, tx_work_dir):
    utils.safe_makedir(tx_work_dir)
    cmd = [sys.executable, os.path.realpath(utils.which("configureStrelkaSomaticWorkflow.py"))]
    cmd += ["--referenceFasta=%s" % ref_file,
            "--callRegions=%s" % _get_region_bed(region, [paired.tumor_data, paired.normal_data], out_file),
            "--runDir=%s" % tx_work_dir,
            "--normalBam=%s" % paired.normal_bam, "--tumorBam=%s" % paired.tumor_bam]
    if dd.get_coverage_interval(paired.tumor_data) not in ["genome"]:
        cmd += ["--targeted"]
    do.run(cmd, "Configure Strelka2 germline calling: %s" % paired.tumor_name)
    return os.path.join(tx_work_dir, "runWorkflow.py")

def _fix_tumor_normal_name(in_file, paired):
    """Replace generic TUMOR NORMAL names in VCF with sample names.
    """
    out_file = in_file.replace(".vcf.gz", "-namefix.vcf.gz")
    if not utils.file_exists(out_file):
        with file_transaction(paired.tumor_data, out_file) as tx_out_file:
            cmd = ("zcat {in_file} | "
                   r"sed 's/NORMAL\tTUMOR/{paired.normal_name}\t{paired.tumor_name}/' | "
                   "bgzip -c > {tx_out_file}")
            do.run(cmd.format(**locals()), "Fix tumor normal names in Strelka2 output")
    return vcfutils.bgzip_and_index(out_file, paired.tumor_data["config"])

def _run_somatic(paired, ref_file, assoc_files, region, out_file):
    if not utils.file_exists(out_file):
        work_dir = "%s-work" % utils.splitext_plus(out_file)[0]
        with file_transaction(paired.tumor_data, work_dir) as tx_work_dir:
            workflow_file = _configure_somatic(paired, ref_file, region, out_file, tx_work_dir)
            _run_workflow(paired.tumor_data, workflow_file, tx_work_dir)
        var_dir = os.path.join(work_dir, "results", "variants")
        vcfutils.combine_variant_files([_fix_tumor_normal_name(os.path.join(var_dir, f), paired)
                                        for f in ["somatic.snvs.vcf.gz", "somatic.indels.vcf.gz"]],
                                       out_file, ref_file, paired.tumor_data["config"], region=region)
    return out_file

def _run_workflow(data, workflow_file, work_dir):
    """Run Strelka2 analysis inside prepared workflow directory.
    """
    utils.remove_safe(os.path.join(work_dir, "workspace"))
    cmd = [sys.executable, workflow_file, "-m", "local", "-j", dd.get_num_cores(data), "--quiet"]
    do.run(cmd, "Run Strelka2: %s" % dd.get_sample_name(data))
    utils.remove_safe(os.path.join(work_dir, "workspace"))
