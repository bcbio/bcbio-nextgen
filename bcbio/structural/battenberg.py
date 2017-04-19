"""Copy number variant calling with cgpBattenberg from Sanger

https://github.com/cancerit/cgpBattenberg
"""
import os

import toolz as tz

from bcbio import utils
from bcbio.bam import ref
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import population, vcfutils

def run(items, background=None):
    """Detect copy number variations from tumor/normal samples using Battenberg.
    """
    paired = vcfutils.get_paired_bams([x["align_bam"] for x in items], items)
    if not paired or not paired.normal_bam:
        logger.warn("Battenberg only works on paired tumor/normal inputs, skipping %s"
                    % dd.get_sample_name(items[0]))
        batout = None
    elif not tz.get_in(["genome_resources", "aliases", "human"], paired.tumor_data):
        logger.warn("Battenberg only works on human data, skipping %s"
                    % dd.get_sample_name(items[0]))
        batout = None
    else:
        batout = _do_run(paired)
        batout["variantcaller"] = "battenberg"
    out = []
    for data in items:
        if batout and dd.get_sample_name(data) == paired.tumor_name:
            if "sv" not in data:
                data["sv"] = []
            data["sv"].append(batout)
        out.append(data)
    return out

def _do_run(paired):
    """Perform Battenberg caling with the paired dataset.

    This purposely does not use a temporary directory for the output
    since Battenberg does smart restarts.
    """
    work_dir = _sv_workdir(paired.tumor_data)
    out = _get_battenberg_out(paired, work_dir)
    ignore_file = os.path.join(work_dir, "ignore_chromosomes.txt")
    if len(_missing_files(out)) > 0:
        ref_file = dd.get_ref_file(paired.tumor_data)
        bat_datadir = os.path.normpath(os.path.join(os.path.dirname(ref_file), os.pardir, "battenberg"))
        ignore_file, gl_file = _make_ignore_file(work_dir, ref_file, ignore_file,
                                                 os.path.join(bat_datadir, "impute", "impute_info.txt"))
        tumor_bam = paired.tumor_bam
        normal_bam = paired.normal_bam
        platform = dd.get_platform(paired.tumor_data)
        genome_build = paired.tumor_data["genome_build"]
        # scale cores to avoid over-using memory during imputation
        cores = max(1, int(dd.get_num_cores(paired.tumor_data) * 0.5))
        gender = {"male": "XY", "female": "XX", "unknown": "L"}.get(population.get_gender(paired.tumor_data))
        if gender == "L":
            gender_str = "-ge %s -gl %s" % (gender, gl_file)
        else:
            gender_str = "-ge %s" % (gender)
        r_export_cmd = utils.get_R_exports()
        local_sitelib = utils.R_sitelib()
        cmd = ("export R_LIBS_USER={local_sitelib} && {r_export_cmd} && "
               "battenberg.pl -t {cores} -o {work_dir} -r {ref_file}.fai "
               "-tb {tumor_bam} -nb {normal_bam} -e {bat_datadir}/impute/impute_info.txt "
               "-u {bat_datadir}/1000genomesloci -c {bat_datadir}/probloci.txt "
               "-ig {ignore_file} {gender_str} "
               "-assembly {genome_build} -species Human -platform {platform}")
        do.run(cmd.format(**locals()), "Battenberg CNV calling")
    assert len(_missing_files(out)) == 0, "Missing Battenberg output: %s" % _missing_files(out)
    out["plot"] = _get_battenberg_out_plots(paired, work_dir)
    out["ignore"] = ignore_file
    return out

def _make_ignore_file(work_dir, ref_file, ignore_file, impute_file):
    """Create input files with chromosomes to ignore and gender loci.
    """
    gl_file = os.path.join(work_dir, "gender_loci.txt")
    chroms = set([])
    with open(impute_file) as in_handle:
        for line in in_handle:
            chrom = line.split()[0]
            chroms.add(chrom)
            if not chrom.startswith("chr"):
                chroms.add("chr%s" % chrom)
    with open(ignore_file, "w") as out_handle:
        for contig in ref.file_contigs(ref_file):
            if contig.name not in chroms:
                out_handle.write("%s\n" % contig.name)
    with open(gl_file, "w") as out_handle:
        for contig in ref.file_contigs(ref_file):
            if contig.name in ["Y", "chrY"]:
                # From https://github.com/cancerit/cgpBattenberg/blob/dev/perl/share/gender/GRCh37d5_Y.loci
                positions = [2934912, 4546684, 4549638, 4550107]
                for pos in positions:
                    out_handle.write("%s\t%s\n" % (contig.name, pos))
    return ignore_file, gl_file

def _missing_files(out):
    missing_files = []
    for key, fname in out.items():
        if not os.path.exists(fname):
            missing_files.append((key, fname))
    return missing_files

def _get_battenberg_out(paired, work_dir):
    out = {}
    for key, ext in [("vrn_file", "battenberg_cn.vcf.gz"),
                     ("subclones", "subclones.txt"),
                     ("contamination", "normal_contamination.txt")]:
        out[key] = os.path.join(work_dir, "%s_%s" % (paired.tumor_name, ext))
    return out

def _get_battenberg_out_plots(paired, work_dir):
    out = {}
    for key, ext in [("copynumberprofile", "copynumberprofile"),
                     ("sunrise_plot", "distance"),
                     ("freq_plot", "Tumor")]:
        plot_file = os.path.join(work_dir, "%s_%s.png" % (paired.tumor_name, ext))
        if os.path.exists(plot_file):
            out[key] = plot_file
    return out

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "battenberg"))
